"""Pure-Python HEALPix module.

This is a Python port from the Typescript code here:
https://github.com/michitaro/healpix
"""
import math

# *** Some generic spherical math helpers

PI2 = 2 * math.pi
PI = math.pi
PI_2 = math.pi / 2
PI_4 = math.pi / 4
PI_8 = math.pi / 8

# Just used as type annotation
XYZ = float, float, float
ZA = float, float

def vec2za(x: float, y: float, z: float) -> ZA:
    r2 = x * x + y * y
    if r2 == 0:
        return -1 if z < 0 else 1, 0
    else:
        a = (math.atan2(y, x) + PI2) % PI2
        z /= math.sqrt(z * z + r2)
        return z, a

def za2vec(z: float, a: float) -> XYZ:
    sin_theta = math.sqrt(1 - z * z)
    X = sin_theta * math.cos(a)
    Y = sin_theta * math.sin(a)
    return X, Y, z

def ang2vec(theta: float, phi: float):
    z = math.cos(theta)
    return za2vec(z, phi)

def vec2ang(v: XYZ):
    z, a = vec2za(v[0], v[1], v[2])
    return math.acos(z), a

def angle(a: XYZ, b: XYZ):
    return 2 * math.asin(math.sqrt(distance2(a, b)) / 2)

def distance2(a: XYZ, b: XYZ):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz

# *** Here we start to implement HEALPix

def order2nside(order: int):
    return 1 << order

def nside2npix(nside: int):
    return 12 * nside * nside

def vec2pix_nest(nside: int, v: XYZ):
    z, a = vec2za(v[0], v[1], v[2])
    return za2pix_nest(nside, z, a)

def vec2pix_ring(nside: int, v: XYZ):
    z, a = vec2za(v[0], v[1], v[2])
    return nest2ring(nside, za2pix_nest(nside, z, a))

def ang2pix_nest(nside: int, theta: float, phi: float):
    z = math.cos(theta)
    return za2pix_nest(nside, z, phi)

def ang2pix_ring(nside: int, theta: float, phi: float):
    z = math.cos(theta)
    return nest2ring(nside, za2pix_nest(nside, z, phi))

def nest2ring(nside: int, ipix: int):
    f, x, y = nest2fxy(nside, ipix)
    return fxy2ring(nside, f, x, y)

def ring2nest(nside: int, ipix: int):
    if nside == 1:
        return ipix

    f, x, y = ring2fxy(nside, ipix)
    return fxy2nest(nside, f, x, y)

def ring2fxy(nside: int, ipix: int):
    polar_lim = 2 * nside * (nside - 1)

    if ipix < polar_lim:  # north polar cap
        i = math.floor((math.sqrt(1 + 2 * ipix) + 1) / 2)
        j = ipix - 2 * i * (i - 1)
        f = math.floor(j / i)
        k = j % i
        x = nside - i + k
        y = nside - 1 - k
        return f, x, y
    elif ipix < polar_lim + 8 * nside * nside: # equatorial belt
        k = ipix - polar_lim
        ring = 4 * nside
        i = nside - math.floor(k / ring)
        s = 1 if i % 2 == 0 else 0
        j = 2 * (k % ring) + s
        jj = j - 4 * nside
        ii = i + 5 * nside - 1
        pp = (ii + jj) / 2
        qq = (ii - jj) / 2
        PP = math.floor(pp / nside)
        QQ = math.floor(qq / nside)
        V = 5 - (PP + QQ)
        H = PP - QQ + 4
        f = 4 * V + (H >> 1) % 4
        x = pp % nside
        y = qq % nside
        return f, x, y
    else: # south polar cap
        p = 12 * nside * nside - ipix - 1
        i = math.floor((math.sqrt(1 + 2 * p) + 1) / 2)
        j = p - 2 * i * (i - 1)
        f = 11 - math.floor(j / i)
        k = j % i
        x = i - k - 1
        y = k
        return f, x, y

def pix2vec_nest(nside: int, ipix: int):
    f, x, y = nest2fxy(nside, ipix)
    t, u = fxy2tu(nside, f, x, y)
    z, a = tu2za(t, u)
    return za2vec(z, a)


def pix2ang_nest(nside: int, ipix: int):
    f, x, y = nest2fxy(nside, ipix)
    t, u = fxy2tu(nside, f, x, y)
    z, a = tu2za(t, u)
    return math.acos(z), a

def pix2vec_ring(nside: int, ipix: int):
    return pix2vec_nest(nside, ring2nest(nside, ipix))

def pix2ang_ring(nside: int, ipix: int):
    return pix2ang_nest(nside, ring2nest(nside, ipix))

def max_pixrad(nside: int):
    unit = PI_4 / nside
    return angle(
        tu2vec(unit, nside * unit),
        tu2vec(unit, (nside + 1) * unit),
    )

def tu2vec(t: float, u: float) -> XYZ:
    z, a = tu2za(t, u)
    return za2vec(z, a)

def corners_nest(nside: int, ipix: int):
    f, x, y = nest2fxy(nside, ipix)
    t, u = fxy2tu(nside, f, x, y)
    d = PI_4 / nside
    xyzs = []
    for (tt, uu) in [
        [0, d],
        [-d, 0],
        [0, -d],
        [d, 0],
    ]:
        z, a = tu2za(t + tt, u + uu)
        xyzs.append(za2vec(z, a))

    return xyzs


def corners_ring(nside: int, ipix: int):
    return corners_nest(nside, ring2nest(nside, ipix))

def nside2pixarea(nside: int):
    return PI / (3 * nside * nside)

def nside2resol(nside: int):
    return math.sqrt(PI / 3) / nside

def pixcoord2vec_nest(nside: int, ipix: int, ne: float, nw: float):
    f, x, y = nest2fxy(nside, ipix)
    t, u = fxy2tu(nside, f, x, y)
    d = PI_4 / nside
    z, a = tu2za(t + d * (ne - nw), u + d * (ne + nw - 1))
    return za2vec(z, a)

def pixcoord2vec_ring(nside: int, ipix: int, ne: float, nw: float):
    return pixcoord2vec_nest(nside, ring2nest(nside, ipix), ne, nw)

def za2pix_nest(nside: int, z: float, a: float):
    t, u = za2tu(z, a)
    f, x, y = tu2fxy(nside, t, u)
    return fxy2nest(nside, f, x, y)

def tu2fxy(nside: int, t: float, u: float):
    f, p, q = tu2fpq(t, u)
    x = clamp(math.floor(nside * p), 0, nside - 1)
    y = clamp(math.floor(nside * q), 0, nside - 1)
    return f, x, y

def sigma(z: float)-> float:
    if z < 0:
        return -sigma(-z)
    else:
        return 2 - math.sqrt(3 * (1 - z))

def za2tu(z: float, a: float):
    """HEALPix spherical projection."""
    if math.fabs(z) <= 2. / 3.: # equatorial belt
        t = a
        u = 3 * PI_8 * z
        return t, u
    else: # polar caps
        p_t = a % (PI_2)
        sigma_z = sigma(z)
        t = a - (math.fabs(sigma_z) - 1) * (p_t - PI_4)
        u = PI_4 * sigma_z
        return t, u

def tu2za(t: float, u: float):
    """Inverse HEALPix spherical projection."""
    abs_u = math.fabs(u)

    if abs_u >= PI_2: # error
        return sign(u), 0
    elif abs_u <= PI_4: # equatorial belt
        z = 8 / (3 * PI) * u
        a = t
        return z, a
    else: #polar caps
        t_t = t % PI_2
        a = t - (abs_u - PI_4) / (abs_u - PI_2) * (t_t - PI_4)
        z = sign(u) * (1 - 1 / 3 * square(2 - 4 * abs_u / PI))
        return z, a

def tu2fpq(t: float, u: float):
    """
    spherical projection -> f, p, q
    f: base pixel index
    p: coord in north east axis of base pixel
    q: coord in north west axis of base pixel
    """
    t /= PI_4
    u /= PI_4
    t = wrap(t, 8)
    t += -4
    u += 5
    pp = clamp((u + t) / 2, 0, 5)
    PP = math.floor(pp)
    qq = clamp((u - t) / 2, 3 - PP, 6 - PP)
    QQ = math.floor(qq)
    V = 5 - (PP + QQ)
    if V < 0: # clamp
        return 0, 1, 1
    H = PP - QQ + 4
    f = 4 * V + (H >> 1) % 4
    p = pp % 1
    q = qq % 1
    return f, p, q

def fxy2nest(nside: int, f: int, x: float, y: float):
    return f * nside * nside + bit_combine(x, y)

def bit_combine(x: int, y: int):
    assert x == int(x)
    assert y == int(y)
    assert x < (1 << 16)
    assert y < (1 << 15)
    x, y = int(x), int(y)

    return (
        x & 1 | (x & 0x2 | y & 0x1) << 1 | (x & 0x4 | y & 0x2) << 2 |
        (x & 0x8 | y & 0x4) << 3 | (x & 0x10 | y & 0x8) << 4 | (x & 0x20 | y & 0x10) << 5 |
        (x & 0x40 | y & 0x20) << 6 | (x & 0x80 | y & 0x40) << 7 | (x & 0x100 | y & 0x80) << 8 |
        (x & 0x200 | y & 0x100) << 9 | (x & 0x400 | y & 0x200) << 10 | (x & 0x800 | y & 0x400) << 11 |
        (x & 0x1000 | y & 0x800) << 12 | (x & 0x2000 | y & 0x1000) << 13 | (x & 0x4000 | y & 0x2000) << 14 |
        (x & 0x8000 | y & 0x4000) << 15 | y & 0x8000 << 16
    )

def bit_decombine(p: int):
    assert p == int(p)
    assert p <= 0x7fffffff
    p = int(p)

    x = (
        (p & 0x1) >> 0 | (p & 0x4) >> 1 | (p & 0x10) >> 2 |
        (p & 0x40) >> 3 | (p & 0x100) >> 4 | (p & 0x400) >> 5 |
        (p & 0x1000) >> 6 | (p & 0x4000) >> 7 | (p & 0x10000) >> 8 |
        (p & 0x40000) >> 9 | (p & 0x100000) >> 10 | (p & 0x400000) >> 11 |
        (p & 0x1000000) >> 12 | (p & 0x4000000) >> 13 |
        (p & 0x10000000) >> 14 | (p & 0x40000000) >> 15
    )
    
    y = (
        (p & 0x2) >> 1 | (p & 0x8) >> 2 | (p & 0x20) >> 3 |
        (p & 0x80) >> 4 | (p & 0x200) >> 5 | (p & 0x800) >> 6 |
        (p & 0x2000) >> 7 | (p & 0x8000) >> 8 | (p & 0x20000) >> 9 |
        (p & 0x80000) >> 10 | (p & 0x200000) >> 11 | (p & 0x800000) >> 12 |
        (p & 0x2000000) >> 13 | (p & 0x8000000) >> 14 | (p & 0x20000000) >> 15
    )

    return x, y


def nest2fxy(nside: int, ipix: int):
    nside2 = nside * nside
    f = math.floor(ipix / nside2)
    k = ipix % nside2 # nested pixel index in base pixel
    x, y = bit_decombine(k)
    return f, x, y

def fxy2ring(nside: int, f: int, x: float, y: float):
    f_row = math.floor(f / 4) # {0 .. 2}
    f1 = f_row + 2            # {2 .. 4}
    v = x + y
    i = f1 * nside - v - 1

    if i < nside: # north polar cap
        f_col = f % 4
        ipix = 2 * i * (i - 1) + (i * f_col) + nside - y - 1
        return ipix
    elif i < 3 * nside: # equatorial belt
        h = x - y
        f2 = 2 * (f % 4) - (f_row % 2) + 1  # {0 .. 7}
        k = (f2 * nside + h + (8 * nside)) % (8 * nside)
        offset = 2 * nside * (nside - 1)
        ipix = offset + (i - nside) * 4 * nside + (k >> 1)
        return ipix
    else: # south polar cap
        i_i = 4 * nside - i
        i_f_col = 3 - (f % 4)
        j = 4 * i_i - (i_i * i_f_col) - y
        i_j = 4 * i_i - j + 1
        ipix = 12 * nside * nside - 2 * i_i * (i_i - 1) - i_j
        return ipix

def fxy2tu(nside: int, f: int, x: float, y: float):
    f_row = math.floor(f / 4)
    f1 = f_row + 2
    f2 = 2 * (f % 4) - (f_row % 2) + 1
    v = x + y
    h = x - y
    i = f1 * nside - v - 1
    k = (f2 * nside + h + (8 * nside))
    t = k / nside * PI_4
    u = PI_2 - i / nside * PI_4
    return t, u

def encode_id(order: int, index: int)-> int:
    return 4 * ((1 << (2 * order)) - 1) + index

def decode_id(id: int) -> (int, int):
    assert id <= 0x7fffffff
    order = 0
    l = (id >> 2) + 1
    while l >= 4:
        l >>= 2
        order += 1

    index = id - (((1 << (2 * order)) - 1) << 2)
    return order, index

def sign(a):
    if a > 0:
        return 1
    elif a < 0:
        return -1
    else:
        return 0

def wrap(a: float, b: float):
    return b - (-a % b) if a < 0 else a % b

def square(a: float):
    return a * a

def clamp(z: float, a: float, b: float):
    return a if z < a else b if z > b else z

# TODO: add query_disc_inclusive_nest
# TODO: add query_disc_inclusive_ring
# TODO: add walk_ring_around
# TODO: add center_t
# TODO: add walk_ring
# TODO: add right_next_pixel

if __name__ == '__main__':
    theta, phi = pix2ang_ring(16, 1440)
    print(theta, PI2 - phi)

    # https://github.com/astropy/astropy-healpix/pull/84
    # This implementation seems to have the same issue
    nside = 2 ** 27
    frac = 2 / 3
    ipix = int(frac * 12 * nside ** 2)
    theta, phi = pix2ang_nest(nside, ipix)
    print('nside:', nside)
    print('ipix=', ipix)
    print('type(theta)=', type(theta))
    print('theta=', theta)
    print('phi=', phi)
    print('distance_to_south_pole=', math.pi - theta)
