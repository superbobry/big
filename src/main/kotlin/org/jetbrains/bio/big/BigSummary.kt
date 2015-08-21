package org.jetbrains.bio.big

import com.google.common.primitives.Doubles
import com.google.common.primitives.Floats
import com.google.common.primitives.Ints
import com.google.common.primitives.Longs

public data class BigSummary(
        /** An upper bound on the number of (bases) with actual data. */
        public var count: Long = 0L,
        /** Minimum item value. */
        public var minValue: Double = java.lang.Double.POSITIVE_INFINITY,
        /** Maximum item value. */
        public var maxValue: Double = java.lang.Double.NEGATIVE_INFINITY,
        /** Sum of values for each base. */
        public var sum: Double = 0.0,
        /** Sum of squares for each base. */
        public var sumSquares: Double = 0.0) {

    /** Returns `true` if a summary contains no data. */
    fun isEmpty() = count == 0L

    fun update(value: Double, intersection: Int, total: Int) {
        val weight = intersection.toDouble() / total
        count += intersection
        sum += value * weight;
        sumSquares += value * value * weight
        minValue = Math.min(minValue, value);
        maxValue = Math.max(maxValue, value);
    }

    /** Because monoids rock. */
    fun plus(other: BigSummary): BigSummary = when {
        isEmpty()       -> other
        other.isEmpty() -> this
        else -> BigSummary(count + other.count,
                           Math.min(minValue, other.minValue),
                           Math.max(maxValue, other.maxValue),
                           sum + other.sum,
                           sumSquares + other.sumSquares)
    }

    fun write(output: OrderedDataOutput) = with(output) {
        writeLong(count)
        writeDouble(minValue)
        writeDouble(maxValue)
        writeDouble(sum)
        writeDouble(sumSquares)
    }

    companion object {
        val BYTES = Longs.BYTES + Doubles.BYTES * 4

        fun read(input: OrderedDataInput) = with(input) {
            val count = readLong()
            val minValue = readDouble()
            val maxValue = readDouble()
            val sum = readDouble()
            val sumSquares = readDouble()
            BigSummary(count, minValue, maxValue, sum, sumSquares)
        }
    }
}

data class ZoomLevel(public val reduction: Int,
                     public val dataOffset: Long,
                     public val indexOffset: Long) {

    fun write(output: OrderedDataOutput) = with(output) {
        writeInt(reduction)
        writeInt(0)  // reserved.
        writeLong(dataOffset)
        writeLong(indexOffset)
    }

    companion object {
        val BYTES = Ints.BYTES * 2 + Longs.BYTES * 2

        fun read(input: OrderedDataInput) = with(input) {
            val reduction = readInt()
            val reserved = readInt()
            check(reserved == 0)
            val dataOffset = readLong()
            val indexOffset = readLong()
            ZoomLevel(reduction, dataOffset, indexOffset)
        }
    }
}

fun List<ZoomLevel>.pick(desiredReduction: Int): ZoomLevel? {
    require(desiredReduction >= 0, "desired must be >=0")
    return if (desiredReduction <= 1) {
        null
    } else {
        var acc = Int.MAX_VALUE
        var closest: ZoomLevel? = null
        for (zoomLevel in this) {
            val d = desiredReduction - zoomLevel.reduction
            if (d >= 0 && d < Math.min(desiredReduction, acc)) {
                acc = d
                closest = zoomLevel
            }
        }

        closest
    }
}

data class ZoomData(
        /** Chromosome id as defined by B+ tree. */
        val chromIx: Int,
        /** 0-based start offset (inclusive). */
        val startOffset: Int,
        /** 0-based end offset (exclusive). */
        val endOffset: Int,
        /**
         * These are just inlined fields of [BigSummary] downcasted
         * to 4 bytes. Top-notch academic design! */
        val count: Int,
        val minValue: Float,
        val maxValue: Float,
        val sum: Float,
        val sumSquares: Float) {

    val interval: ChromosomeInterval get() = Interval(chromIx, startOffset, endOffset)

    /** Returns `true` if zoom contains no data. */
    fun isEmpty() = count == 0

    fun write(output: OrderedDataOutput) = with(output) {
        writeInt(chromIx)
        writeInt(startOffset)
        writeInt(endOffset)
        writeInt(count)
        writeFloat(minValue)
        writeFloat(maxValue)
        writeFloat(sum)
        writeFloat(sumSquares)
    }

    companion object {
        val SIZE: Int = Ints.BYTES * 3 + Ints.BYTES + Floats.BYTES * 4

        fun read(input: OrderedDataInput): ZoomData = with(input) {
            val chromIx = readInt()
            val startOffset = readInt()
            val endOffset = readInt()
            val count = readInt()
            val minValue = readFloat()
            val maxValue = readFloat()
            val sum = readFloat();
            val sumSquares = readFloat();
            return ZoomData(chromIx, startOffset, endOffset, count,
                            minValue, maxValue, sum, sumSquares);
        }
    }
}

fun Pair<ChromosomeInterval, BigSummary>.toZoomData(): ZoomData {
    val (chromIx, startOffset, endOffset) = first
    val (count, minValue, maxValue, sum, sumSquares) = second
    return ZoomData(chromIx, startOffset, endOffset,
                    count.toInt(),
                    minValue.toFloat(), maxValue.toFloat(),
                    sum.toFloat(), sumSquares.toFloat())
}