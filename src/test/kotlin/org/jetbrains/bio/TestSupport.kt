package org.jetbrains.bio

import org.apache.commons.lang3.SystemUtils
import org.jetbrains.bio.big.BigBedFile
import org.jetbrains.bio.big.BigFile
import org.jetbrains.bio.big.BigWigFile
import java.io.IOException
import java.nio.ByteOrder
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

internal object Examples {
    @JvmStatic operator fun get(name: String): Path {
        val url = Examples.javaClass.classLoader.getResource(name)
                  ?: throw IllegalStateException("resource not found")

        return Paths.get(url.toURI()).toFile().toPath()
    }
}

/** Fetches chromosome sizes from a UCSC provided TSV file. */
internal fun Path.chromosomes(): List<Pair<String, Int>> {
    return bufferedReader().lineSequence().map { line ->
        val chunks = line.split('\t', limit = 3)
        chunks[0] to chunks[1].toInt()
    }.toList()
}

internal inline fun withTempFile(prefix: String, suffix: String,
                                 block: (Path) -> Unit) {
    val path = Files.createTempFile(prefix, suffix)
    try {
        block(path)
    } finally {
        try {
            Files.delete(path)
        } catch (e: IOException) {
            // Mmaped buffer not yet garbage collected. Leave it to the VM.
            path.toFile().deleteOnExit()
        }
    }
}

abstract class NamedRomBufferFactoryProvider(private val title: String) {
    abstract operator fun invoke(path: String, byteOrder: ByteOrder, limit: Long = -1L): RomBufferFactory
    override fun toString() = title
}

private val mmapSupported =
    SystemUtils.OS_ARCH in arrayOf("amd64", "x86_64") && (SystemUtils.IS_OS_LINUX || SystemUtils.IS_OS_MAC_OSX)

fun romFactoryProviders(): List<NamedRomBufferFactoryProvider> {
    val providers: MutableList<NamedRomBufferFactoryProvider> = mutableListOf(
            object : NamedRomBufferFactoryProvider("EndianSynchronizedBufferFactory") {
                override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) =
                        EndianSynchronizedBufferFactory.create(path, byteOrder)
            },

            object : NamedRomBufferFactoryProvider("EndianBufferFactory") {
                override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) =
                        EndianBufferFactory.create(path, byteOrder)
            },

            object : NamedRomBufferFactoryProvider("EndianThreadSafeBufferFactory") {
                override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) =
                        EndianThreadSafeBufferFactory(path, byteOrder)
            }
    )
    if (mmapSupported) {
        providers.add(object : NamedRomBufferFactoryProvider("MMBRomBufferFactory") {
            override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) = MMBRomBufferFactory(Paths.get(path), byteOrder)
        })
    }
    return providers
}

fun threadSafeRomFactoryProvidersAndPrefetchParams(): List<Array<Any>> {
    val providers: MutableList<NamedRomBufferFactoryProvider> = mutableListOf(
            object : NamedRomBufferFactoryProvider("EndianThreadSafeBufferFactory") {
                override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) =
                        EndianThreadSafeBufferFactory(path, byteOrder)
            },
            object : NamedRomBufferFactoryProvider("EndianSynchronizedBufferFactory") {
                override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) =
                        EndianSynchronizedBufferFactory.create(path, byteOrder)
            }
    )
    if (mmapSupported) {
        providers.add(object : NamedRomBufferFactoryProvider("MMBRomBufferFactory") {
            override fun invoke(path: String, byteOrder: ByteOrder, limit: Long) = MMBRomBufferFactory(Paths.get(path), byteOrder)
        })
    }
    return providers.flatMap { fp ->
        intArrayOf(BigFile.PREFETCH_LEVEL_OFF, BigFile.PREFETCH_LEVEL_FAST, BigFile.PREFETCH_LEVEL_DETAILED).map { prefetch ->
            arrayOf(fp, prefetch)
        }
    }
}

fun romFactoryProviderParams(): List<Array<Any>> = romFactoryProviders().map { arrayOf<Any>(it) }

fun romFactoryProviderAndPrefetchParams(): List<Array<Any>> = romFactoryProviders().flatMap { fp ->
    intArrayOf(BigFile.PREFETCH_LEVEL_OFF, BigFile.PREFETCH_LEVEL_FAST, BigFile.PREFETCH_LEVEL_DETAILED).map { prefetch ->
        arrayOf(fp, prefetch)
    }
}

fun BigBedFile.Companion.read(path: Path, provider: NamedRomBufferFactoryProvider,
                              prefetch: Int) =
        BigBedFile.read(path.toString(), prefetch) { src, byteOrder ->
            provider(src, byteOrder)
        }

fun BigWigFile.Companion.read(path: Path, provider: NamedRomBufferFactoryProvider,
                              prefetch: Int) =
        BigWigFile.read(path.toString(), prefetch) { src, byteOrder ->
            provider(src, byteOrder)
        }

fun BigFile.Companion.read(path: Path, provider: NamedRomBufferFactoryProvider,
                           prefetch: Int) =
        BigFile.read(path.toString(), prefetch) { src, byteOrder ->
            provider(src, byteOrder)
        }

fun romFactoryByteOrderCompressionParamsSets(): Iterable<Array<Any>> {
    val params = mutableListOf<Array<Any>>()
    for (factoryProvider in romFactoryProviders()) {
        for (byteOrder in listOf(ByteOrder.BIG_ENDIAN, ByteOrder.LITTLE_ENDIAN)) {
            CompressionType.values().mapTo(params) { arrayOf(byteOrder, it, factoryProvider) }
        }
    }
    return params
}

fun allBigFileParams(): List<Array<Any>> = romFactoryByteOrderCompressionParamsSets().flatMap {
    intArrayOf(BigFile.PREFETCH_LEVEL_OFF, BigFile.PREFETCH_LEVEL_FAST, BigFile.PREFETCH_LEVEL_DETAILED).map { prefetch ->
        it + prefetch
    }
}
