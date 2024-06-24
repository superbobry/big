package org.jetbrains.bio.big

import com.google.common.collect.ComparisonChain
import com.google.common.io.Closeables
import org.jetbrains.bio.*
import java.io.IOException
import java.nio.ByteOrder
import java.nio.file.Path
import java.util.*

/**
 * Just like BED only BIGGER.
 */
class BigBedFile private constructor(
        path: String,
        buffFactory: RomBufferFactory,
        magic: Int,
        prefetch: Int,
        cancelledChecker: (() -> Unit)?
) : BigFile<BedEntry>(path, buffFactory, magic, prefetch, cancelledChecker) {

    override fun summarizeInternal(
            input: RomBuffer,
            query: ChromosomeInterval,
            numBins: Int,
            cancelledChecker: (() -> Unit)?
    ): Sequence<IndexedValue<BigSummary>> {
        val coverage = query(input, query, overlaps = true, cancelledChecker = cancelledChecker).aggregate()
        var edge = 0
        return query.slice(numBins).mapIndexed { i, bin ->
            val summary = BigSummary()
            for (j in edge until coverage.size) {
                val (bedEntry, score) = coverage[j]
                if (bedEntry.end <= bin.startOffset) {
                    edge = j + 1
                    continue
                } else if (bedEntry.start > bin.endOffset) {
                    break
                }

                val interval = Interval(query.chromIx, bedEntry.start, bedEntry.end)
                if (interval intersects bin) {
                    summary.update(score.toDouble(),
                            interval.intersectionLength(bin))
                }
            }

            if (summary.isEmpty()) null else IndexedValue(i, summary)
        }.filterNotNull()
    }

    /**
     * Returns `true` if a given entry is consistent with the query.
     * That is
     *   it either intersects the query (and overlaps is `true`)
     *   or it is completely contained in the query.
     */
    private fun ChromosomeInterval.contains(startOffset: Int, endOffset: Int,
                                            overlaps: Boolean): Boolean {
        val interval = Interval(chromIx, startOffset, endOffset)
        return (overlaps && interval intersects this) || interval in this
    }

    override fun queryInternal(decompressedBlock: RomBuffer,
                               query: ChromosomeInterval,
                               overlaps: Boolean): Sequence<BedEntry> {
        val chrom = chromosomes[query.chromIx]

        return with(decompressedBlock) {
            val chunk = ArrayList<BedEntry>()
            do {
                val chromIx = readInt()
                assert(chromIx == query.chromIx) {
                    "interval contains wrong chromosome $chromIx, expected ${query.chromIx}, source: $source"
                }
                val startOffset = readInt()
                val endOffset = readInt()
                val rest = readCString()
                if (query.contains(startOffset, endOffset, overlaps)) {
                    chunk.add(BedEntry(chrom, startOffset, endOffset, rest))
                }
            } while (hasRemaining())

            chunk.asSequence()
        }
    }

    companion object {
        /** Magic number used for determining [ByteOrder]. */
        internal const val MAGIC = 0x8789F2EB.toInt()

        @Throws(IOException::class)
        @JvmStatic
        fun read(path: Path, cancelledChecker: (() -> Unit)? = null) =
                read(path.toString(), cancelledChecker = cancelledChecker)

        @Throws(IOException::class)
        @JvmStatic
        fun read(src: String, prefetch: Int = PREFETCH_LEVEL_DETAILED,
                 cancelledChecker: (() -> Unit)? = null,
                 factoryProvider: RomBufferFactoryProvider = defaultFactory()
        ): BigBedFile {
            val factory = factoryProvider(src, ByteOrder.LITTLE_ENDIAN)
            try {
                val byteOrder = getByteOrder(src, MAGIC, factory)
                factory.order = byteOrder

                return BigBedFile(src, factory, MAGIC, prefetch, cancelledChecker)
            } catch (e: Exception) {
                Closeables.close(factory, true)
                throw e
            }
        }

        private class BedEntrySummary {
            val chromosomes = HashSet<String>()
            var count = 0
            var sum = 0L

            /** Makes sure the entries are sorted by offset. */
            private var edge = 0

            /** Makes sure the entries are sorted by chromosome. */
            private var previous = ""

            operator fun invoke(entry: BedEntry) {
                val switch = entry.chrom !in chromosomes
                require(entry.chrom == previous || switch) {
                    "must be sorted by chromosome"
                }

                require(entry.start >= edge || switch) { "must be sorted by offset" }

                chromosomes.add(entry.chrom)
                sum += entry.end - entry.start
                count++

                previous = entry.chrom
                edge = entry.start
            }
        }

        /**
         * Creates a BigBED file from given entries.
         *
         * @param bedEntries entries sorted by chromosome *and* start offset.
         *                   The method traverses the entries twice:
         *                   firstly to summarize and secondly to write
         *                   and index.
         * @param chromSizes chromosome names and sizes, e.g.
         *                   `("chrX", 59373566)`. Entries on chromosomes
         *                   missing from this list will be dropped.
         * @param outputPath BigBED file path.
         * @param itemsPerSlot number of items to store in a single
         *                     R+ tree index node. Defaults to `1024`.
         * @param zoomLevelCount number of zoom levels to pre-compute.
         *                       Defaults to `8`.
         * @param compression method for data sections, see [CompressionType].
         * @param order byte order used, see [java.nio.ByteOrder].
         * @param cancelledChecker Throw cancelled exception to abort operation
         * @@throws IOException if any of the read or write operations failed.
         */
        @Throws(IOException::class)
        @JvmStatic
        @JvmOverloads
        fun write(
                bedEntries: Iterable<BedEntry>,
                chromSizes: Iterable<Pair<String, Int>>,
                outputPath: Path,
                itemsPerSlot: Int = 1024, zoomLevelCount: Int = 8,
                compression: CompressionType = CompressionType.SNAPPY,
                order: ByteOrder = ByteOrder.nativeOrder(),
                cancelledChecker: (() -> Unit)? = null) {

            val summary = BedEntrySummary().apply { bedEntries.forEach { this(it) } }

            val header = OrderedDataOutput(outputPath, order).use { output ->
                output.skipBytes(Header.BYTES)
                output.skipBytes(ZoomLevel.BYTES * zoomLevelCount)
                val totalSummaryOffset = output.tell()
                output.skipBytes(BigSummary.BYTES)

                val unsortedChromosomes = chromSizes.filter { it.first in summary.chromosomes }
                        .mapIndexed { i, (key, size) -> BPlusLeaf(key, i, size) }
                val chromTreeOffset = output.tell()
                BPlusTree.write(output, unsortedChromosomes)

                val unzoomedDataOffset = output.tell()
                val resolver = unsortedChromosomes.map { it.key to it.id }.toMap()
                val leaves = ArrayList<RTreeIndexLeaf>()
                var uncompressBufSize = 0

                cancelledChecker?.invoke()
                for ((chrName, items) in bedEntries.asSequence().groupingBy { it.chrom }) {
                    cancelledChecker?.invoke()

                    val chromIx = resolver[chrName]
                    if (chromIx == null) {
                        items.forEach {}  // Consume.
                        continue
                    }

                    val it = items.iterator()

                    while (it.hasNext()) {
                        val dataOffset = output.tell()
                        var leafStart = Int.MAX_VALUE
                        var leafEnd = 0
                        val current = output.with(compression) {
                            for ((_, start, end, rest) in it.asSequence().take(itemsPerSlot)) {
                                writeInt(chromIx)
                                writeInt(start)
                                writeInt(end)
                                writeString(rest)
                                writeByte(0)  // NUL-terminated.

                                leafStart = Math.min(leafStart, start)
                                leafEnd = Math.max(leafEnd, end)
                            }
                        }

                        leaves.add(RTreeIndexLeaf(
                                Interval(chromIx, leafStart, leafEnd),
                                dataOffset, output.tell() - dataOffset))
                        uncompressBufSize = Math.max(uncompressBufSize, current)
                    }
                }

                cancelledChecker?.invoke()
                val unzoomedIndexOffset = output.tell()
                RTreeIndex.write(output, leaves, itemsPerSlot = itemsPerSlot)

                Header(
                        output.order, MAGIC,
                        version = if (compression == CompressionType.SNAPPY) 5 else 4,
                        zoomLevelCount = zoomLevelCount,
                        chromTreeOffset = chromTreeOffset,
                        unzoomedDataOffset = unzoomedDataOffset,
                        unzoomedIndexOffset = unzoomedIndexOffset,
                        fieldCount = 3, definedFieldCount = 3,
                        totalSummaryOffset = totalSummaryOffset,
                        uncompressBufSize = if (compression.absent) 0 else uncompressBufSize)
            }

            OrderedDataOutput(outputPath, order, create = false).use { header.write(it) }

            with(summary) {
                if (count > 0) {
                    val initial = Math.max(sum divCeiling count.toLong(), 1).toInt() * 10
                    Post.zoom(outputPath, itemsPerSlot, initial = initial, cancelledChecker = cancelledChecker)
                }
            }

            Post.totalSummary(outputPath)
        }
    }
}

private class AggregationEvent(val offset: Int, val type: Int,
                               val item: BedEntry) : Comparable<AggregationEvent> {

    override fun toString() = "${if (type == END) "END" else "START"}@$offset"

    override fun compareTo(other: AggregationEvent): Int = ComparisonChain.start()
            .compare(offset, other.offset)
            .compare(type, other.type)
            .result()
}

private const val END = 0    // must be before start.
private const val START = 1

/** Computes intervals of uniform coverage. */
internal fun Sequence<BedEntry>.aggregate(): List<Pair<BedEntry, Int>> {
    val events = flatMap {
        sequenceOf(AggregationEvent(it.start, START, it),
                AggregationEvent(it.end, END, it))
    }.toMutableList()

    events.sort()

    var current = 0
    var left = 0
    val res = ArrayList<Pair<BedEntry, Int>>()
    for ((i, event) in events.withIndex()) {
        when {
            event.type == START -> {
                if (current == 0) {
                    left = event.offset
                }

                current += 1
            }
            event.type == END || i == events.size - 1 -> {
                assert(event.offset >= left)
                // Produce a single aggregate for duplicate intervals.
                // For ease of use we abuse the semantics of the
                // '#score' field in 'BedEntry'.
                // This method is required for coverage summarising, so
                // additional fields after '#score' not used, let's skip them
                if (event.offset > left) {
                    val item = event.item
                    res.add(BedEntry(item.chrom, left, event.offset, "") to current)

                    left = event.offset
                }

                current -= 1
            }
        }
    }

    return res
}
