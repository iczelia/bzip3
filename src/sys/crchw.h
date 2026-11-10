/*  Copyright (C) 2022-2026 Kamila Szewczyk

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
    General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/*  The body of a CRC-32C kernel built on a CPU instruction.  The includer defines crc_w, CRC_W, step(),
    CRC_BYTE and CRC_FN.  Three independent lanes keep the instruction, which has a latency of three cycles,
    busy every cycle; the lanes are joined by advancing a checksum over the bytes of the lanes after it.  */

static INLINE crc_w lanes(crc_w c, const u8 ** pp, sz * pn, sz len,
                          const u32 z[4][256]) {
  const u8 * p = *pp;
  sz n = *pn;
  for (; n >= 3 * len; n -= 3 * len, p += 2 * len) {
    crc_w c1 = 0, c2 = 0;
    const u8 * e = p + len;
    do {
      c = step(c, p);
      c1 = step(c1, p + len);
      c2 = step(c2, p + 2 * len);
      p += CRC_W;
    } while (p < e);
    c = bz3_crc_shift(z, bz3_crc_shift(z, (u32) c) ^ (u32) c1) ^ (u32) c2;
  }
  *pp = p;  *pn = n;
  return c;
}

u32 CRC_FN(const bz3_crc * t, u32 crc, const void * buf, sz n) {
  const u8 * p = buf;
  crc_w c = crc;
  for (; n && ((uintptr_t) p & (CRC_W - 1)); n--)
    c = CRC_BYTE((u32) c, *p++);
  c = lanes(c, &p, &n, CRC_LONG, t->u.z.zlong);
  c = lanes(c, &p, &n, CRC_SHORT, t->u.z.zshort);
  for (; n >= CRC_W; n -= CRC_W, p += CRC_W) c = step(c, p);
  for (; n; n--) c = CRC_BYTE((u32) c, *p++);
  return (u32) c;
}
