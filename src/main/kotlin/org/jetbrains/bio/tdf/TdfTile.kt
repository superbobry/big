package org.jetbrains.bio.tdf

import org.apache.log4j.Logger
import org.jetbrains.bio.RomBuffer
import org.jetbrains.bio.ScoredInterval

/**
 * Data container in [TdfFile].
 *
 * @since 0.2.2
 */
interface TdfTile {
    fun view(trackNumber: Int) = TdfTileView(this, trackNumber)

    fun getValue(trackNumber: Int, idx: Int): Float

    fun getStartOffset(idx: Int): Int

    fun getEndOffset(idx: Int): Int

    /** Number of data points in a tile. */
    val size: Int

    companion object {
        private val LOG = Logger.getLogger(TdfTile::class.java)

        internal fun read(input: RomBuffer, expectedTracks: Int) = with(input) {
            val type = getCString()
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
                val start = tile.getStartOffset(it)
                val end = tile.getEndOffset(it)
                ScoredInterval(start, end, value)
            }
        }.iterator()
    }
}

data class TdfBedTile(val starts: IntArray, val ends: IntArray,
                      val data: Array<FloatArray>) : TdfTile {
    override val size: Int get() = starts.size

    override fun getStartOffset(idx: Int) = starts[idx]

    override fun getEndOffset(idx: Int) = ends[idx]

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: RomBuffer, expectedTracks: Int) = with(input) {
            val size = getInt()

            val start = getInts(size)
            val end = getInts(size)

            val trackCount = getInt()
            check(trackCount == expectedTracks) {
                "expected $expectedTracks tracks, got: $trackCount"
            }

            val data = Array(trackCount) {
                getFloats(size)
            }

            TdfBedTile(start, end, data)
        }
    }
}

data class TdfFixedTile(val start: Int, val span: Double,
                        val data: Array<FloatArray>) : TdfTile {

    override val size: Int get() = data.first().size

    override fun getStartOffset(idx: Int): Int {
        return start + (idx * span).toInt()
    }

    override fun getEndOffset(idx: Int): Int {
        return start + ((idx + 1) * span).toInt()
    }

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: RomBuffer, expectedTracks: Int) = with(input) {
            val size = getInt()
            val start = getInt()
            val span = getInt().toDouble()

            // vvv not part of the implementation, see igvteam/igv/#180.
            // val trackCount = readInt()
            val data = Array(expectedTracks) {
                getFloats(size)
            }

            TdfFixedTile(start, span, data)
        }
    }
}

data class TdfVaryTile(val starts: IntArray,
                       val span: Int,
                       val data: Array<FloatArray>) : TdfTile {

    override val size: Int get() = starts.size

    override fun getStartOffset(idx: Int) = starts[idx]

    override fun getEndOffset(idx: Int) = (starts[idx] + span).toInt()

    override fun getValue(trackNumber: Int, idx: Int) = data[trackNumber][idx]

    companion object {
        fun fill(input: RomBuffer, expectedTracks: Int) = with(input) {
            // This is called 'tiledStart' in IGV sources and is unused.
            getInt() // start
            val span = getFloat().toInt()  // Really?
            val size = getInt()

            val step = getInts(size)

            val trackCount = getInt()
            check(trackCount == expectedTracks) {
                "expected $expectedTracks tracks, got: $trackCount"
            }
            
            val data = Array(trackCount) {
                getFloats(size)
            }

            TdfVaryTile(step, span, data)
        }
    }
}