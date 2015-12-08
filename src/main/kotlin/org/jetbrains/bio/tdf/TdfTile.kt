package org.jetbrains.bio.tdf

import org.apache.log4j.Logger
import org.jetbrains.bio.OrderedDataInput
import org.jetbrains.bio.ScoredInterval

/**
 * Data container in [TdfFile].
 *
 * @since 0.2.2
 */
interface TdfTile {
    fun view(trackNumber: Int) = TdfTileView(this, trackNumber)

    fun getValue(trackNumber: Int, idx: Int): Float

    // TODO: position implies 1-based indexing, which isn't the case.
    fun getStartPosition(idx: Int): Int

    fun getEndPosition(idx: Int): Int

    /** Number of data points in a tile. */
    val size: Int

    companion object {
        private val LOG = Logger.getLogger(TdfTile::class.java)

        internal fun read(input: OrderedDataInput, expectedTracks: Int) = with(input) {
            val type = readCString()
            when (type) {
                "fixedStep" -> TdfFixedTile.fill(this, expectedTracks)
                "variableStep" -> TdfVaryTile.fill(this, expectedTracks)
                "bed", "bedWithName" -> {
                    if (type === "bedWithName") {
                        LOG.warn("bedWithName is not supported, assuming bed")
                    }

                    TdfBedTile.fill(this, expectedTracks)
                }
                else -> error("unexpected type: $type")
            }
        }
    }
}

/** A view of a single track in a [TdfTile]. */
class TdfTileView(private val tile: TdfTile,
                  private val trackNumber: Int) : Iterable<ScoredInterval?> {
    override fun iterator(): Iterator<ScoredInterval?> {
        return (0 until tile.size).asSequence().map {
            val value = tile.getValue(trackNumber, it)
            if (value.isNaN()) {
                null
            } else {
                val start = tile.getStartPosition(it)
                val end = tile.getEndPosition(it)
                ScoredInterval(start, end, value)
            }
        }.iterator()
    }
}

data class TdfBedTile(val starts: IntArray, val ends: IntArray,
                      val data: Array<FloatArray>) : TdfTile {
    override val size: Int get() = starts.size

    override fun getStartPosition(idx: Int) = starts[idx]

    override fun getEndPosition(idx: Int) = ends[idx]

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: OrderedDataInput, expectedTracks: Int) = with(input) {
            val size = readInt()
            val start = IntArray(size)
            for (i in 0 until size) {
                start[i] = readInt()
            }
            val end = IntArray(size)
            for (i in 0 until size) {
                end[i] = readInt()
            }

            val trackCount = readInt()
            check(trackCount == expectedTracks) {
                "expected $expectedTracks tracks, got: $trackCount"
            }
            val data = Array(trackCount) {
                val acc = FloatArray(size)
                for (i in 0 until size) {
                    acc[i] = readFloat()
                }
                acc
            }

            TdfBedTile(start, end, data)
        }
    }
}

data class TdfFixedTile(val start: Int, val span: Double,
                        val data: Array<FloatArray>) : TdfTile {

    override val size: Int get() = data.first().size

    override fun getStartPosition(idx: Int): Int {
        return start + (idx * span).toInt()
    }

    override fun getEndPosition(idx: Int): Int {
        return start + ((idx + 1) * span).toInt()
    }

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: OrderedDataInput, expectedTracks: Int) = with(input) {
            val size = readInt()
            val start = readInt()
            val span = readFloat().toDouble()

            // vvv not part of the implementation, see igvteam/igv/#180.
            // val trackCount = readInt()
            val data = Array(expectedTracks) {
                val acc = FloatArray(size)
                for (i in 0 until size) {
                    acc[i] = readFloat()
                }
                acc
            }

            TdfFixedTile(start, span, data)
        }
    }
}

data class TdfVaryTile(val starts: IntArray, val span: Int,
                       val data: Array<FloatArray>) : TdfTile {

    override val size: Int get() = starts.size

    override fun getStartPosition(idx: Int) = starts[idx]

    override fun getEndPosition(idx: Int) = (starts[idx] + span).toInt()

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: OrderedDataInput, expectedTracks: Int) = with(input) {
            // This is called 'tiledStart' in IGV sources and is unused.
            val start = readInt()
            val span = readFloat().toInt()  // Really?
            val size = readInt()

            val step = IntArray(size)
            for (i in 0 until size) {
                step[i] = readInt()
            }

            val trackCount = readInt()
            check(trackCount == expectedTracks) {
                "expected $expectedTracks tracks, got: $trackCount"
            }
            val data = Array(trackCount) {
                val acc = FloatArray(size)
                for (i in 0 until size) {
                    acc[i] = readFloat()
                }

                acc
            }

            TdfVaryTile(step, span, data)
        }
    }
}