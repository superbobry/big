package org.jbb.big;

import java.io.IOException;
import java.nio.ByteOrder;

/**
 * Chromosome B+ tree headers for bigBED.
 *
 * @author Sergey Zherevchuk
 */
public class BptHeader {

  public static final int MAGIC = 0x78ca8c91;

  public ByteOrder byteOrder = ByteOrder.BIG_ENDIAN;
  public final int blockSize;
  public final int keySize;
  public final int valSize;
  public final long itemCount;
  public final long rootOffset;

  public BptHeader(final ByteOrder byteOrder, final int blockSize, final  int keySize,
                   final int valSize, final long itemCount, final long rootOffset) {
    this.byteOrder = byteOrder;
    this.blockSize = blockSize;
    this.keySize = keySize;
    this.valSize = valSize;
    this.itemCount = itemCount;
    this.rootOffset = rootOffset;
  }

  // Attach B+ tree of chromosome names and ids.
  public static BptHeader read(final SeekableStream s, final long chromTreeOffset)
      throws IOException {
    s.seek(chromTreeOffset);
    s.guess(MAGIC);
    final int blockSize = s.readInt();
    final int keySize = s.readInt();
    final int valSize = s.readInt();
    final long itemCount = s.readLong();
    s.readLong(); // reserved bytes
    final long rootOffset = s.tell();

    return new BptHeader(s.order(), blockSize, keySize, valSize, itemCount, rootOffset);
  }
}
