# =================================================================
# 1. Global Cimport
# =================================================================

from libc.stdint cimport uint8_t, uint64_t
from libc.string cimport memcpy, memset

# =================================================================

# =================================================================
# 2. CONST UINT64
# =================================================================

cdef uint64_t PI    = 0x243F6A8885A308D3
cdef uint64_t PHI   = 0x9E3779B97F4A7C15
cdef uint64_t ZETA  = 0x73356E847241A417
cdef uint64_t THETA = 0x486445D7E84B8E1D
cdef uint64_t E     = 0xB7E151628AED2A6A
cdef uint64_t SQRT2 = 0x6A09E667F3BCC908
cdef uint64_t LN2   = 0x317208E206385C69
cdef uint64_t CBRT3 = 0x502859495143A675
cdef uint64_t SQRT7 = 0x51E2838C61E1A334
cdef uint64_t SQRT8 = 0x6A09E667F3BCC909
cdef uint64_t CBRT7 = 0x6338B57065997232
cdef uint64_t CBRT9 = 0x787D6E90F845C22D

cdef uint64_t Bit0Mask = 0xFFFFFFFFFFFFFFFE
# =================================================================

# =================================================================
# 3. GLOBAL FUNCTION
# =================================================================

cdef inline uint64_t mix(const uint64_t i, const uint64_t c, const uint64_t idx):
    cdef uint64_t x = i + idx
    x ^= x >> 32
    x *= idx
    x ^= x >> 32
    x *= c
    x ^= x >> 32
    return x

cdef inline uint64_t rol64(const uint64_t x, const uint64_t c):
    cdef uint64_t n = c & 63
    return (x << n) | (x >> (64 - n))

cdef inline uint64_t ror64(const uint64_t x, const uint64_t c):
    cdef uint64_t n = c & 63
    return (x >> n) | (x << (64 - n))

# =================================================================

# =================================================================
# 4. DRA2
# =================================================================

# =================================================================
# 4.1 DRA2-128
# =================================================================

cdef class DRA2_128:
    cdef uint64_t[2] state
    cdef uint64_t[16] buffer
    cdef uint64_t[2] _digest
    cdef size_t length
    cdef size_t total

    def __cinit__(self):
        self.state[0] = SQRT8
        self.state[1] = CBRT9
        self.length = 0
        self.total = 0
    
    cdef void absorb(self):
        cdef uint64_t a = self.state[0], b = self.state[1]
        for r in range(0, 8, 1):
            a ^= self.buffer[r]
            a += PI
            a *= PI
        for r in range(8, 16, 1):
            b ^= self.buffer[r]
            b += PHI
            b *= PHI
        self.state[0] = a
        self.state[1] = b
        return

    cdef void permute(self):
        cdef uint64_t a = self.state[0], b = self.state[1]

        for _ in range(5):
            a ^= PI
            b ^= PHI
            a *= rol64(b * mix(b, CBRT3, 128), 11)
            b *= rol64(a * mix(a, E, 256)    , 13)

        cdef uint64_t final_a = (a >> 33) * (b >> 33)
        cdef uint64_t final_b = (b >> 33) * (a >> 33) 

        self._digest[0] = final_a
        self._digest[1] = final_b

    cpdef DRA2_128 update(self, const uint8_t[:] obj):
        if obj.shape[0] <= 0:
            return self
        cdef size_t remaining = obj.shape[0]
        cdef const uint8_t *ptr = &obj[0]
        cdef size_t BLOCKSIZE = 128
        cdef size_t to_cpy
        self.total += remaining

        while remaining > 0:
            to_cpy = min(remaining, BLOCKSIZE - self.length)
            memcpy(<uint8_t*>self.buffer + self.length, ptr, to_cpy)
            self.length += to_cpy
            ptr += to_cpy
            remaining -= to_cpy
            
            if self.length == BLOCKSIZE:
                self.absorb()
                self.length = 0
        return self
    
    cpdef DRA2_128 copy(self):
        cdef DRA2_128 newobj = DRA2_128.__new__(DRA2_128)
        newobj.state = self.state
        newobj.buffer = self.buffer
        newobj.length = self.length
        newobj._digest = self._digest
        newobj.total = self.total
        return newobj
    
    cpdef bytes digest(self):
        cdef uint8_t *b_ptr = <uint8_t*>self.buffer
        if self.length > 0:
            memset(b_ptr + self.length, 0x80, 1)
            memset(b_ptr + self.length + 1, 0x01, 120 - self.length)
        (<uint64_t*>(b_ptr + 120))[0] = self.total
        self.absorb()
        
        self.permute()
        return (<uint8_t*>self._digest)[:16]
    
    cpdef str hexdigest(self):
        return self.digest().hex()

cpdef DRA2_128 dra2_128(const uint8_t[:] obj = None):
    cdef DRA2_128 newobj = DRA2_128()
    if obj != None:
        newobj.update(obj)
    
    return newobj

# =================================================================

# =================================================================
# 4.2 DRA2-224
# =================================================================

cdef class DRA2_224:
    cdef uint64_t[4] state
    cdef uint64_t[32] buffer
    cdef uint64_t[4] _digest
    cdef size_t length
    cdef size_t total

    def __cinit__(self):
        self.state[0] = SQRT2
        self.state[1] = SQRT7
        self.state[2] = CBRT3
        self.state[3] = ZETA
        self.length = 0
        self.total = 0
    
    cdef void absorb(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]
        
        for i in range(0, 8, 1):
            a ^= self.buffer[i]
            a *= E

        for i in range(8, 16, 1):
            b ^= self.buffer[i]
            b *= SQRT2
        
        for i in range(16, 24, 1):
            c ^= self.buffer[i]
            c *= LN2
        
        for i in range(24, 32, 1):
            d ^= self.buffer[i]
            d *= CBRT3

        self.state[0] = a
        self.state[1] = b
        self.state[2] = c
        self.state[3] = d
        return
    
    cdef void permute(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]

        for _ in range(5):
            a ^= PI
            b ^= PHI
            c ^= ZETA
            d ^= THETA

            a *= rol64(d * mix(d, CBRT3, 512), 15)
            b *= rol64(c * mix(c, E, 1024)    , 17)
            c *= rol64(b * mix(b, LN2, 2048)  , 19)
            d *= rol64(a * mix(a, THETA, 4096), 21)

        cdef uint64_t final_a = (c >> 33) * (d >> 33)
        cdef uint64_t final_b = (d >> 33) * (c >> 33)
        cdef uint64_t final_c = (a >> 33) * (b >> 33)
        cdef uint64_t final_d = (b >> 33) * (a >> 33)

        self._digest[0] = final_a
        self._digest[1] = final_b
        self._digest[2] = final_c
        self._digest[3] = final_d ^ (final_d >> 32)

    cpdef DRA2_224 update(self, const uint8_t[:] obj):
        if obj.shape[0] <= 0:
            return self
        cdef size_t remaining = obj.shape[0]
        cdef const uint8_t *ptr = &obj[0]
        cdef size_t BLOCKSIZE = 128
        cdef size_t to_cpy
        self.total += remaining

        while remaining > 0:
            to_cpy = min(remaining, BLOCKSIZE - self.length)
            memcpy(<uint8_t*>self.buffer + self.length, ptr, to_cpy)
            self.length += to_cpy
            ptr += to_cpy
            remaining -= to_cpy
            if self.length == BLOCKSIZE:
                self.absorb()
                self.length = 0
        return self

    cpdef DRA2_224 copy(self):
        cdef DRA2_224 newobj = DRA2_224.__new__(DRA2_224)
        newobj.state = self.state
        newobj.buffer = self.buffer
        newobj.length = self.length
        newobj._digest = self._digest
        newobj.total = self.total
        return newobj

    cpdef bytes digest(self):
        cdef uint8_t *b_ptr = <uint8_t*>self.buffer
        if self.length > 0:
            memset(b_ptr + self.length, 0x80, 1)
            memset(b_ptr + self.length + 1, 0x01, 120 - self.length)
        (<uint64_t*>(b_ptr + 120))[0] = self.total
        self.absorb()
        
        self.permute()
        return (<uint8_t*>self._digest)[:28]

    cpdef str hexdigest(self):
        return self.digest().hex()

cpdef DRA2_224 dra2_224(const uint8_t[:] obj = None):
    cdef DRA2_224 newobj = DRA2_224()

    if obj != None:
        newobj.update(obj)
    
    return newobj

# =================================================================

# =================================================================
# 4.3 DRA2-256
# =================================================================
cdef class DRA2_256:
    cdef uint64_t[4] state
    cdef uint64_t[36] buffer
    cdef uint64_t[4] _digest
    cdef size_t length
    cdef size_t total

    def __cinit__(self):
        self.state[0] = PI
        self.state[1] = PHI
        self.state[2] = ZETA
        self.state[3] = THETA
        self.length = 0
        self.total = 0
    
    cdef void absorb(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]
        for i in range(0, 9, 1):
            a ^= self.buffer[i]
            a += PI
            a ^= PI
        for i in range(9, 18, 1):
            b ^= self.buffer[i]
            b += PHI
            b ^= PHI
        for i in range(18, 27, 1):
            c ^= self.buffer[i]
            c += ZETA
            c ^= ZETA
        for i in range(27, 36, 1):
            d ^= self.buffer[i]
            d += THETA
            d ^= THETA
        
        self.state[0] = a
        self.state[1] = b
        self.state[2] = c
        self.state[3] = d
        return
    
    cdef void permute(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]

        for _ in range(5):
            a ^= PI
            b ^= PHI
            c ^= ZETA
            d ^= THETA

            a *= rol64(d * mix(d, CBRT3, 8192), 23)
            b *= rol64(c * mix(c, E, 16384)    , 25)
            c *= rol64(b * mix(b, LN2, 32968)  , 27)
            d *= rol64(a * mix(a, THETA, 65936), 29)

        cdef uint64_t final_a = (b >> 32) * (c >> 32)
        cdef uint64_t final_b = (c >> 32) * (d >> 32)
        cdef uint64_t final_c = (d >> 32) * (a >> 32)
        cdef uint64_t final_d = (a >> 32) * (b >> 32)

        self._digest[0] = final_a
        self._digest[1] = final_b
        self._digest[2] = final_c
        self._digest[3] = final_d
    
    cpdef DRA2_256 update(self, const uint8_t[:] obj):
        if obj.shape[0] <= 0:
            return self
        cdef size_t remaining = obj.shape[0]
        cdef const uint8_t *ptr = &obj[0]
        cdef size_t BLOCKSIZE = 288
        cdef size_t to_cpy
        self.total += remaining

        while remaining > 0:
            to_cpy = min(remaining, BLOCKSIZE - self.length)
            memcpy(<uint8_t*>self.buffer + self.length, ptr, to_cpy)
            self.length += to_cpy
            ptr += to_cpy
            remaining -= to_cpy
            
            if self.length == BLOCKSIZE:
                self.absorb()
                self.length = 0
        return self
    
    cpdef DRA2_256 copy(self):
        cdef DRA2_256 newobj = DRA2_256.__new__(DRA2_256)
        newobj.state = self.state
        newobj.buffer = self.buffer
        newobj.length = self.length
        newobj._digest = self._digest
        newobj.total = self.total
        return newobj

    cpdef bytes digest(self):
        cdef uint8_t *b_ptr = <uint8_t*>self.buffer
        if self.length > 0:
            memset(b_ptr + self.length, 0x80, 1)
            memset(b_ptr + self.length + 1, 0x01, 280 - self.length)
        (<uint64_t*>(b_ptr + 280))[0] = self.total
        self.absorb()
        
        self.permute()
        return (<uint8_t*>self._digest)[:32]
    
    cpdef str hexdigest(self):
        return self.digest().hex()

cpdef DRA2_256 dra2_256(const uint8_t[:] obj = None):
    cdef DRA2_256 newobj = DRA2_256()

    if obj != None:
        newobj.update(obj)
    
    return newobj
# =================================================================

# =================================================================
# 4.4 DRA2-384
# =================================================================


# =================================================================

# =================================================================
# 4.5 DRA2-512
# =================================================================
# =================================================================

# =================================================================
# 4.6 DRAGON2
# =================================================================

cdef class DRAGON2:
    cdef uint64_t[4] state
    cdef uint64_t[40] buffer
    cdef uint64_t[4] _digest
    cdef size_t length
    cdef size_t total

    def __cinit__(self):
        self.state[0] = PI
        self.state[1] = PHI
        self.state[2] = ZETA
        self.state[3] = THETA
        self.length = 0
        self.total = 0
    
    cdef void absorb(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]
        for i in range(0, 10, 1):
            a ^= self.buffer[i]
            a += PI
            
        for i in range(10, 20, 1):
            b ^= self.buffer[i]
            b += PHI
            
        for i in range(20, 30, 1):
            c ^= self.buffer[i]
            c += ZETA
        for i in range(30, 40, 1):
            d ^= self.buffer[i]
            d += THETA
        
        self.state[0] = a
        self.state[1] = b
        self.state[2] = c
        self.state[3] = d
    
    cdef void absorb_direct(self, const uint8_t* ptr):
        # Cast con trỏ về uint64_t để đọc 8-byte một lần (tốc độ tối đa)
        cdef uint64_t* data64 = <uint64_t*>ptr
        cdef size_t i

        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]
        for i in range(0, 10, 1):
            a ^= data64[i]
            a += PI

        for i in range(10, 20, 1):
            b ^= data64[i]
            b += PHI

        for i in range(20, 30, 1):
            c ^= data64[i]
            c += ZETA

        for i in range(30, 40, 1):
            d ^= data64[i]
            d += THETA

        
        self.state[0] = a
        self.state[1] = b
        self.state[2] = c
        self.state[3] = d
    cdef void permute(self):
        cdef uint64_t a = self.state[0], b = self.state[1], c = self.state[2], d = self.state[3]

        for i in range(5):
            a ^= PI
            b ^= PHI
            c ^= ZETA
            d ^= THETA

            a *= rol64(a * mix(d, CBRT3, THETA), a)
            b *= rol64(b * mix(c, E, LN2)      , b)
            c *= rol64(c * mix(b, LN2, E)      , c)
            d *= rol64(d * mix(a, THETA, CBRT3), d)

        cdef uint64_t final_a = (c >> 33) * (d >> 33)
        cdef uint64_t final_b = (d >> 33) * (c >> 33)
        cdef uint64_t final_c = (a >> 33) * (b >> 33)
        cdef uint64_t final_d = (b >> 33) * (a >> 33)

        self._digest[0] = final_a
        self._digest[1] = final_b
        self._digest[2] = final_c
        self._digest[3] = final_d

    cpdef DRAGON2 update(self, const uint8_t[:] obj):
        cdef size_t remaining = obj.shape[0]
        if remaining <= 0:
            return self
        cdef const uint8_t *ptr = &obj[0]
        cdef size_t BLOCKSIZE = 320
        cdef size_t to_cpy
        self.total += remaining

        while remaining > 0:
            to_cpy = min(remaining, BLOCKSIZE - self.length)
            memcpy(<uint8_t*>self.buffer + self.length, ptr, to_cpy)
            self.length += to_cpy
            ptr += to_cpy
            remaining -= to_cpy
            if self.length == BLOCKSIZE:
                self.absorb()
                self.length = 0
        return self
    cpdef DRAGON2 copy(self):
        cdef DRAGON2 newobj = DRAGON2.__new__(DRAGON2)
        newobj.state = self.state
        newobj.buffer = self.buffer
        newobj.length = self.length
        newobj._digest = self._digest
        newobj.total = self.total
        return newobj

    cpdef bytes digest(self, size_t digest_size):
        cdef size_t outputsize = min(32, max(8, digest_size))

        cdef uint8_t *b_ptr = <uint8_t*>self.buffer
        if self.length > 0:
            memset(b_ptr + self.length, 0x80, 1)
            memset(b_ptr + self.length + 1, 0x01, 312 - self.length)
        (<uint64_t*>(b_ptr + 312))[0] = self.total ^ outputsize
        self.absorb()
        
        self.permute()
        
        return (<uint8_t*>self._digest)[:outputsize]

    cpdef str hexdigest(self, size_t digest_size):
        return self.digest(digest_size).hex()

cpdef DRAGON2 dragon2(const uint8_t[:] obj = None):
    cdef DRAGON2 newobj = DRAGON2()
    if obj != None:
        newobj.update(obj)
    return newobj

# =================================================================