package org.jetbrains.bio.big

import com.google.common.collect.Lists
import com.google.common.primitives.Ints
import java.io.IOException
import java.nio.file.Path
import java.util.Collections
import kotlin.platform.platformStatic

/**
 * Bigger brother of the good-old WIG format.
 */
public class BigWigFile throws(IOException::class) protected constructor(path: Path) :
        BigFile<WigSection>(path, magic = BigWigFile.MAGIC) {

    override fun summarizeInternal(query: ChromosomeInterval, numBins: Int): List<BigSummary> {
        val wigItems = query(query).flatMap { it.query().asSequence() }.toList()
        var edge = 0
        return query.slice(numBins).map { bin ->
            val summary = BigSummary()
            for (j in edge until wigItems.size()) {
                val wigItem = wigItems[j]
                if (wigItem.end <= bin.startOffset) {
                    edge = j + 1
                    continue
                } else if (wigItem.start > bin.endOffset) {
                    break
                }

                val interval = Interval(query.chromIx, wigItem.start, wigItem.end)
                if (interval intersects bin) {
                    summary.update(wigItem.score.toDouble(),
                                   (interval intersection bin).length(),
                                   interval.length())
                }
            }

            summary
        }.toList()
    }

    override fun queryInternal(dataOffset: Long, dataSize: Long,
                               query: ChromosomeInterval): Sequence<WigSection> {
        val chrom = chromosomes[query.chromIx]
        return sequenceOf(input.with(dataOffset, dataSize, compressed) {
            val chromIx = readInt()
            assert(chromIx == query.chromIx, "section contains wrong chromosome")
            val start = readInt()
            readInt()   // end.
            val step = readInt()
            val span = readInt()
            val type = readUnsignedByte()
            readByte()  // reserved.
            val count = readUnsignedShort()

            val types = WigSection.Type.values()
            check(type >= 1 && type <= types.size())
            when (types[type - 1]) {
                WigSection.Type.BED_GRAPH ->
                    throw IllegalStateException("bedGraph sections aren't supported")
                WigSection.Type.VARIABLE_STEP -> {
                    val section = VariableStepSection(chrom, span)
                    for (i in 0 until count) {
                        val position = readInt()
                        section[position] = readFloat()
                    }

                    section
                }
                WigSection.Type.FIXED_STEP -> {
                    val section = FixedStepSection(chrom, start, step, span)
                    for (i in 0 until count) {
                        section.add(readFloat())
                    }

                    section
                }
            }
        })
    }

    companion object {
        /** Magic number used for determining [ByteOrder]. */
        private val MAGIC: Int = 0x888FFC26.toInt()

        throws(IOException::class)
        public platformStatic fun read(path: Path): BigWigFile = BigWigFile(path)

        /**
         * Creates a BigWIG file from given sections.
         *
         * @param wigSections sections to write and index.
         * @param chromSizesPath path to the TSV file with chromosome
         *                       names and sizes.
         * @param outputPath BigWIG file path.
         * @param compressed compress BigWIG data sections with gzip.
         *                   Defaults to `false`.
         * @throws IOException if any of the read or write operations failed.
         */
        throws(IOException::class)
        public platformStatic fun write(wigSections: Iterable<WigSection>,
                                        chromSizesPath: Path,
                                        outputPath: Path,
                                        compressed: Boolean = true) {
            SeekableDataOutput.of(outputPath).use { output ->
                output.skipBytes(0, BigFile.Header.BYTES)

                val unsortedChromosomes = chromSizesPath.bufferedReader()
                        .lineSequence().mapIndexed { i, line ->
                    val chunks = line.split('\t', limit = 3)
                    BPlusLeaf(chunks[0], i, chunks[1].toInt())
                }.toList()

                val chromTreeOffset = output.tell()
                BPlusTree.write(output, unsortedChromosomes)

                // XXX move to 'BedFile'?
                val unzoomedDataOffset = output.tell()
                val resolver = unsortedChromosomes.map { it.key to it.id }.toMap()
                val leaves = Lists.newArrayList<RTreeIndexLeaf>()
                var uncompressBufSize = 0
                wigSections.groupBy { it.chrom }.forEach { entry ->
                    val (name, sections) = entry
                    Collections.sort(sections) { e1, e2 -> Ints.compare(e1.start, e2.start) }

                    val chromId = resolver[name]!!
                    for (section in sections.asSequence().flatMap { it.splice() }) {
                        val dataOffset = output.tell()
                        val current = output.with(compressed) {
                            when (section) {
                                is FixedStepSection -> section.write(this, resolver)
                                is VariableStepSection -> section.write(this, resolver)
                            }
                        }

                        leaves.add(RTreeIndexLeaf(Interval(chromId, section.start, section.end),
                                                  dataOffset, output.tell() - dataOffset))
                        uncompressBufSize = Math.max(uncompressBufSize, current)
                    }
                }

                val unzoomedIndexOffset = output.tell()
                RTreeIndex.write(output, leaves, itemsPerSlot = 1)

                val header = BigFile.Header(
                        output.order,
                        chromTreeOffset = chromTreeOffset,
                        unzoomedDataOffset = unzoomedDataOffset,
                        unzoomedIndexOffset = unzoomedIndexOffset,
                        fieldCount = 0, definedFieldCount = 0,
                        uncompressBufSize = if (compressed) uncompressBufSize else 0)
                header.write(output, MAGIC)
            }
        }
    }
}

private fun FixedStepSection.write(output: OrderedDataOutput, resolver: Map<String, Int>) {
    with(output) {
        writeInt(resolver[chrom]!!)
        writeInt(start)
        writeInt(end)
        writeInt(step)
        writeInt(span)
        writeByte(WigSection.Type.FIXED_STEP.ordinal() + 1)
        writeByte(0) // reserved.
        writeShort(values.size())
        for (i in 0 until values.size()) {
            writeFloat(values[i])
        }
    }
}

private fun VariableStepSection.write(output: OrderedDataOutput, resolver: Map<String, Int>) {
    with(output) {
        writeInt(resolver[chrom]!!)
        writeInt(start)
        writeInt(end)
        writeInt(0)
        writeInt(span)
        writeByte(WigSection.Type.VARIABLE_STEP.ordinal() + 1)
        writeByte(0)  // reserved.
        writeShort(values.size())
        for (i in 0 until values.size()) {
            writeInt(positions[i])
            writeFloat(values[i])
        }
    }
}