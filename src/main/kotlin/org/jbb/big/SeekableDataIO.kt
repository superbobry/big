package org.jbb.big

import com.google.common.primitives.Ints
import com.google.common.primitives.Longs
import com.google.common.primitives.Shorts
import java.io.*
import java.nio.ByteOrder
import java.nio.channels.Channels
import java.nio.file.Path
import java.util.zip.Deflater
import java.util.zip.DeflaterOutputStream
import java.util.zip.Inflater
import java.util.zip.InflaterInputStream

/**
 * A byte order-aware seekable complement to [java.io.DataInputStream]
 * and [java.io.DataOutputStream].
 *
 * @author Sergei Lebedev
 * @since 29/06/15
 */
public open class SeekableDataInput(private val path: Path,
                                    public var order: ByteOrder) :
        DataInput, Closeable, AutoCloseable {
    protected val file: RandomAccessFile = RandomAccessFile(path.toFile(), "r")

    protected open val input: DataInput get() = file

    /** Guess byte order from a given big-endian `magic`. */
    public fun guess(magic: Int) {
        val b = ByteArray(4)
        readFully(b)
        val bigMagic = Ints.fromBytes(b[0], b[1], b[2], b[3])
        order = if (bigMagic != magic) {
            val littleMagic = Ints.fromBytes(b[3], b[2], b[1], b[0])
            check(littleMagic == magic, "bad signature")
            ByteOrder.LITTLE_ENDIAN
        } else {
            ByteOrder.BIG_ENDIAN
        }
    }

    public fun compressed(size: Long): CompressedDataInput {
        return CompressedDataInput(path, tell(), order, size)
    }

    override fun readFully(b: ByteArray?) = input.readFully(b)

    override fun readFully(b: ByteArray?, off: Int, len: Int) = input.readFully(b, off, len)

    override fun skipBytes(n: Int): Int = input.skipBytes(n)

    public fun seek(pos: Long): Unit = file.seek(pos)

    public fun tell(): Long = file.getFilePointer()

    override fun readBoolean(): Boolean = input.readBoolean()

    override fun readByte(): Byte = input.readByte()

    override fun readUnsignedByte(): Int = input.readUnsignedByte()

    override fun readChar(): Char = input.readChar()

    override fun readShort(): Short {
        val b1 = readAndCheckByte()
        val b2 = readAndCheckByte()
        return if (order == ByteOrder.BIG_ENDIAN) {
            Shorts.fromBytes(b1, b2)
        } else {
            Shorts.fromBytes(b2, b1)
        }
    }

    override fun readUnsignedShort(): Int {
        val b1 = readAndCheckByte()
        val b2 = readAndCheckByte()
        return if (order == ByteOrder.BIG_ENDIAN) {
            Ints.fromBytes(0, 0, b1, b2)
        } else {
            Ints.fromBytes(0, 0, b2, b1)
        }
    }

    override fun readInt(): Int {
        val b1 = readAndCheckByte()
        val b2 = readAndCheckByte()
        val b3 = readAndCheckByte()
        val b4 = readAndCheckByte()
        return if (order == ByteOrder.BIG_ENDIAN) {
            Ints.fromBytes(b1, b2, b3, b4)
        } else {
            Ints.fromBytes(b4, b3, b2, b1)
        }
    }

    override fun readLong(): Long {
        val b1 = readAndCheckByte()
        val b2 = readAndCheckByte()
        val b3 = readAndCheckByte()
        val b4 = readAndCheckByte()
        val b5 = readAndCheckByte()
        val b6 = readAndCheckByte()
        val b7 = readAndCheckByte()
        val b8 = readAndCheckByte()
        return if (order == ByteOrder.BIG_ENDIAN) {
            Longs.fromBytes(b1, b2, b3, b4, b5, b6, b7, b8)
        } else {
            Longs.fromBytes(b8, b7, b6, b5, b4, b3, b2, b1)
        }
    }

    override fun readFloat(): Float = java.lang.Float.intBitsToFloat(readInt())

    override fun readDouble(): Double = java.lang.Double.longBitsToDouble(readLong());

    override fun readLine(): String? = throw UnsupportedOperationException()

    override fun readUTF(): String = throw UnsupportedOperationException()

    throws(IOException::class)
    private fun readAndCheckByte(): Byte {
        val b = file.read()
        if (b == -1) {
            throw EOFException()
        }

        return b.toByte()
    }

    override fun close() = file.close()

    companion object {
        public fun of(path: Path,
                      order: ByteOrder = ByteOrder.BIG_ENDIAN): SeekableDataInput {
            return SeekableDataInput(path, order)
        }
    }
}

public open class SeekableDataOutput(private val path: Path,
                                     public var order: ByteOrder) :
        DataOutput, Closeable, AutoCloseable {
    protected val file: RandomAccessFile = RandomAccessFile(path.toFile(), "rw")

    protected open val output: DataOutput get() = file

    public fun compressed(size: Long): CompressedDataOutput {
        return CompressedDataOutput(path, tell(), order, size)
    }

    public fun seek(pos: Long): Unit = file.seek(pos)

    public fun tell(): Long = file.getFilePointer()

    override fun write(b: ByteArray?) = output.write(b)

    override fun write(b: ByteArray?, off: Int, len: Int) = output.write(b, off, len)

    override fun writeBytes(s: String) = output.writeBytes(s)

    public fun writeBytes(s: String, length: Int) {
        file.writeBytes(s)
        writeByte(0, length - s.length())
    }

    override fun writeBoolean(v: Boolean) = output.writeBoolean(v)

    override fun writeByte(v: Int) = output.writeByte(v)

    public fun writeByte(v: Int, count: Int) {
        for (i in 0 until count) {
            writeByte(v)
        }
    }

    override fun write(b: Int) = output.write(b)

    override fun writeChar(v: Int) = output.writeChar(v)

    override fun writeShort(v: Int) {
        output.writeShort(if (order == ByteOrder.BIG_ENDIAN) {
            v
        } else {
            val b = Shorts.toByteArray(v.toShort())
            Shorts.fromBytes(b[1], b[0]).toInt()
        })
    }

    public fun writeUnsignedShort(v: Int) {
        output.writeShort(if (order == ByteOrder.BIG_ENDIAN) {
            v
        } else {
            val b = Ints.toByteArray(v)
            Shorts.fromBytes(b[3], b[2]).toInt()
        })
    }

    override fun writeInt(v: Int) {
        output.writeInt(if (order == ByteOrder.BIG_ENDIAN) {
            v
        } else {
            val b = Ints.toByteArray(v)
            Ints.fromBytes(b[3], b[2], b[1], b[0])
        })
    }

    override fun writeLong(v: Long): Unit {
        output.writeLong(if (order == ByteOrder.BIG_ENDIAN) {
            v
        } else {
            val b = Longs.toByteArray(v)
            Longs.fromBytes(b[7], b[6], b[5], b[4], b[3], b[2], b[1], b[0])
        })
    }

    override fun writeFloat(v: Float) = writeInt(java.lang.Float.floatToIntBits(v))

    override fun writeDouble(v: Double): Unit = writeLong(java.lang.Double.doubleToLongBits(v))

    override fun writeChars(s: String): Unit = throw UnsupportedOperationException()

    override fun writeUTF(s: String): Unit = throw UnsupportedOperationException()

    override fun close() = file.close()

    companion object {
        throws(FileNotFoundException::class)
        public fun of(path: Path, order: ByteOrder = ByteOrder.BIG_ENDIAN): SeekableDataOutput {
            return SeekableDataOutput(path, order)
        }
    }
}

class CompressedDataInput(path: Path, offset: Long, order: ByteOrder,
                          private val size: Long) : SeekableDataInput(path, order) {
    private val inflater = Inflater()

    init {
        file.seek(offset)
    }

    override val input: DataInputStream = DataInputStream(
            InflaterInputStream(Channels.newInputStream(file.getChannel()),
                                inflater))

    override fun close() {
        check(inflater.getBytesRead() == size)
        input.close()  // should close the file.
    }
}

class CompressedDataOutput(path: Path, offset: Long, order: ByteOrder,
                           private val size: Long) :
        SeekableDataOutput(path, order) {
    private val deflater = Deflater()

    init {
        file.seek(offset)
    }

    override val output: DataOutputStream = DataOutputStream(
            DeflaterOutputStream(Channels.newOutputStream(file.getChannel()),
                                 deflater))

    override fun close() {
        check(deflater.getBytesWritten() == size)
        output.close()  // should close the file.
    }
}