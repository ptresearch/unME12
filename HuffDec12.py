import os, sys, struct, zlib

class Error(Exception): pass

def cwDec(w): # Convert 16-bit value to string codeword
  return bin(0x10000 | w).rstrip('0')[3:-1]

def cwEnc(cw): # Convert string codeword to 16-bit value
  return int((cw+'1').ljust(16, '0'), 2)

#***************************************************************************
#***************************************************************************
#***************************************************************************

def HuffTabReader_bin(ab):
  fmtRec = struct.Struct("<HB")
  o = 0
  while o < len(ab):
    w, cb = fmtRec.unpack_from(ab, o)
    o += fmtRec.size
    v = ab[o:o+cb]
    assert len(v) == cb
    o += cb
    yield(cwDec(w), cb, v)

def HuffTabPack_bin(d):
  r = []
  for cw in sorted(d.keys())[::-1]:
    v = d[cw]
    if v is None: continue
    r.append(struct.pack("<HB", cwEnc(cw), len(v)) + v)
  return "".join(r)

#***************************************************************************
#***************************************************************************
#***************************************************************************

def HuffTabReader_text(ab):
  for ln in ab.splitlines():
    a = ln.strip().split()
    if len(a) < 1: continue # Skip empty lines
    cw = a[0] # String representation of Huffman codeword
    v = a[1] if len(a) > 1 else None # Huffman sequence value
    if v is None: # Not defined
      cb = None
    elif v.startswith("??"): # Sequence length is known
      cb = len(v)/2
      v = None
    else: # Value is known
      v = v.decode("hex")
      cb = len(v)
    yield(cw, cb, v)

def HuffTabPack_text(dLen, d, mode=None):
  if mode is None: mode = HuffDecoder.DUMP_ALL
  r = []
  for cw in sorted(dLen.keys())[::-1]:
    cb = dLen[cw]
    if cb is None and mode in (HuffDecoder.DUMP_KNOWN, HuffDecoder.DUMP_LEN): continue # Ignore if sequence length is not known
    v = d.get(cw, None)
    if v is None:
      if HuffDecoder.DUMP_KNOWN == mode: continue # Ignore if sequence is not known
      v = "" if cb is None else "??"*cb
    else: v = v.encode("hex").upper()
    pad = '\t'*(2 - len(cw)/8)
    r.append("%s%s%s" % (cw, pad, v))
  return "\n".join(r)

#***************************************************************************
#***************************************************************************
#***************************************************************************

def HuffTab_extendLen(dLen, extShape=False):
  shape = []
  aCW = sorted(dLen.keys())[::-1]
  minBits, maxBits = len(aCW[0]), len(aCW[-1])
  aCW.append('0'*(maxBits+1)) # Longer than max

  nBits = minBits # Current length
  e = int(aCW[0], 2)|1 # End value for current length
  for o in xrange(1, len(aCW)):
    nextBits = len(aCW[o])
    if nextBits == nBits: continue # Run until length change
    assert nextBits > nBits # Length must increase
    s = int(aCW[o-1], 2) # Start value for current length
    for i in xrange(s, e+1):
      cw = bin(i)[2:].zfill(nBits)
      if cw not in dLen: dLen[cw] = None
    e = int(aCW[o], 2)|1 # End value for next length
    shape.append([(x << (16-nBits)) for x in xrange(s, e/2, -1)])
    if extShape:
      for i in xrange(s-1, e/2, -1):
        cw = bin(i)[2:].zfill(nBits)
        if cw not in dLen: dLen[cw] = None
    nBits = nextBits
  return shape

#***************************************************************************
#***************************************************************************
#***************************************************************************

class HuffNode(object):
  def __init__(self, cw, hd):
    self.cw = cw # String codeword value
#    self.w = cwEnc(cw) # Encoded codeword value
    if hd:
      self.nBits = len(cw) # Length of codeword in bits
      self.cb = hd.dLen.get(cw, None)
      self.av = [d.get(cw, None) for d in hd.adTab]
    else:
      self.nBits = None # Actual length of codeword is unknown

#***************************************************************************
#***************************************************************************
#***************************************************************************

class HuffDecoder(object):
  NAMES = ("Code", "Data")
  DUMP_KNOWN = 0
  DUMP_LEN = 1
  DUMP_ALL = 2
  dPrefix = {DUMP_KNOWN:"kno", DUMP_LEN:"len", DUMP_ALL:"all"}
  fmtInt = struct.Struct("<L")
  baseDir = os.path.split(__file__)[0]
  BLOCK_SIZE = 0x1000 # 4K bytes

  def __init__(self, ver=None):
    names = self.NAMES if ver is None else tuple("%s%d" % (n, ver) for n in self.NAMES)
    try:
      self.loadTables(names) # Load from text version
#      with open("huff11.bin", "wb") as fo: fo.write(zlib.compress(self.packTables(), 9)[2:-4])
    except:
      with open(os.path.join(self.baseDir, "huff11.bin"), "rb") as fi: self.unpackTables(zlib.decompress(fi.read(), -15)) # Load from compressed version
    self.prepareMap()

  def loadTable(self, items):
    sv = set() # Set for values
    d = {}
    for cw, cb, v in items:
      if cw in d: raise Error("Codeword %s already defined" % cw)

      if cb is None: continue
      cbKnown = self.dLen.get(cw, None)
      if cbKnown is None: self.dLen[cw] = cb
      elif cb != cbKnown: raise Error("Codeword %s sequence length %d != known %d" % (cw, cb, cbKnown))

      if v is None: continue
      assert len(v) == cb
      d[cw] = v # Remember value
#      if v in sv: raise Error("Value %s already present" % v.encode("hex"))
      sv.add(v)

    self.adTab.append(d)

  def unpackTables(self, ab):
    n, = self.fmtInt.unpack_from(ab)
    o = self.fmtInt.size
    self.dLen, self.adTab = {}, []
    for i in xrange(n):
      cb, = self.fmtInt.unpack_from(ab, o)
      o += self.fmtInt.size
      data = ab[o:o+cb]
      assert len(data) == cb
      o += cb
      self.loadTable(HuffTabReader_bin(data))

  def packTables(self):
    r = [self.fmtInt.pack(len(self.adTab))]
    for d in self.adTab:
      ab = HuffTabPack_bin(d)
      r.append(self.fmtInt.pack(len(ab)) + ab)
    return "".join(r)

  def loadTables(self, names=None):
    if names is None: names = self.NAMES
    self.dLen, self.adTab = {}, []
    for name in names:
      with open(os.path.join(self.baseDir, "%s.txt" % name)) as fi:
        self.loadTable(HuffTabReader_text(fi.read()))

  def saveTables(self, mode=None, names=None):
    if mode is None: mode = self.DUMP_ALL
    if names is None: names = self.NAMES
    dLen = self.dLen.copy()
    HuffTab_extendLen(dLen)
    for name,d in zip(names, self.adTab):
      with open(os.path.join(self.baseDir, "%s%s.txt" % (self.dPrefix[mode], name)), "w") as fo:
        print >>fo, HuffTabPack_text(dLen, d, mode)

  def propagateMap(self, node):
    cw = node.cw
    for idx in xrange(int(cw[::-1], 2), len(self.aMap), 1<<len(cw)):
      assert self.aMap[idx] is None
      self.aMap[idx] = node

  def prepareMap(self):
    aCW = sorted(self.dLen.keys())[::-1]
    minBits, maxBits = len(aCW[0]), len(aCW[-1])
    self.mask = (1 << maxBits) - 1
    self.aMap = [None]*(1<<maxBits) # 2**maxBits map
    aCW.append('0'*(maxBits+1)) # Longer than max
    nBits = minBits # Current length
    e = int(aCW[0], 2)|1 # End value for current length
    for o in xrange(1, len(aCW)):
      nextBits = len(aCW[o])
      if nextBits == nBits: continue # Run until length change
      assert nextBits > nBits # Length must increase
      s = int(aCW[o-1], 2) # Start value for current length
      for i in xrange(s, e+1):
        cw = bin(i)[2:].zfill(nBits)
        self.propagateMap(HuffNode(cw, self))
      e = int(aCW[o], 2)|1 # End value for next length
      for i in xrange(e/2 + 1, s): # Handle values with unknown codeword length
        cw = bin(i)[2:].zfill(nBits)
        self.propagateMap(HuffNode(cw, None))
      nBits = nextBits
    for v in self.aMap: assert v is not None

  def enumCW(self, ab):
    v = int(bin(int("01"+ab.encode("hex"), 16))[3:][::-1], 2) # Reversed bits
    cb = 0
    while cb < self.BLOCK_SIZE: # Block length
      node = self.aMap[v & self.mask]
      if node.nBits is None: raise Error("Unknown codeword %s* length" % node.cw)
      yield node
      v >>= node.nBits
      if node.cb is not None: cb += node.cb

  def decompressChunk(self, ab, iTab):
    r = []
    cb = 0
    for node in self.enumCW(ab):
      v = node.av[iTab]
      if v is None: raise Error("Unknown sequence for codeword %s in table #%d" % (node.cw, iTab))
      r.append(v)
      cb += len(v)
      if cb >= self.BLOCK_SIZE: break
    return "".join(r)

  def decompress(self, ab, length):
    nChunks, left = divmod(length, self.BLOCK_SIZE)
    assert 0 == left
    aOfs = list(struct.unpack_from("<%dL" % nChunks, ab))
    aOpt = [0]*nChunks
    for i in xrange(nChunks):
      aOpt[i], aOfs[i] = divmod(aOfs[i], 0x40000000)

    base = nChunks*4
    aOfs.append(len(ab) - base)
    r = []
    for i, opt in enumerate(aOpt):
      iTab, bCompr = divmod(opt, 2)
      assert 1 == bCompr
      unpacked = self.decompressChunk(ab[base + aOfs[i]: base + aOfs[i+1]], iTab)
      assert len(unpacked) == self.BLOCK_SIZE
      r.append(unpacked)
    return "".join(r)

#***************************************************************************
#***************************************************************************
#***************************************************************************

def main(argv):
  hd = HuffDecoder()
#  with open("huff11.bin", "wb") as fo: fo.write(zlib.compress(hd.packTables(), 9)[2:-4])
#  with open("all.bin", "wb") as fo: fo.write(hd.packTables())
#  for mode in (HuffDecoder.DUMP_KNOWN, HuffDecoder.DUMP_LEN, HuffDecoder.DUMP_ALL): hd.saveTables(mode)

  hd.prepareMap()
  for fn in argv[1:]:
    with open(fn, "rb") as fi: ab = fi.read()
    nChunks, = struct.unpack_from("<L", ab)
    ab = hd.decompress(ab[4:], nChunks * hd.BLOCK_SIZE)
    with open(os.path.splitext(fn)[0] + ".mod", "wb") as fo: fo.write(ab)

if __name__=="__main__": main(sys.argv)
