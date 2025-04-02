import sys
import os
import ctypes

# sage imports
from sage.all import (PolynomialRing)

# FPGA env
path = os.getenv('PATH')
sys.path.append(path)

import utilities.fpga_io as fio

# FPGA driver imports
libdma = ctypes.CDLL("/home/wzy/k-nakamura/richelot_isogeny/driver/tools/libdma.so")

libdma.read_to_buffer.argtypes = [
    ctypes.c_char_p,   # device
    ctypes.c_int,      # fpga_fd
    ctypes.c_char_p,   # buffer
    ctypes.c_uint64,   # size
    ctypes.c_uint64,   # addr
]

libdma.write_from_buffer.argtypes = [
    ctypes.c_char_p,   # device
    ctypes.c_int,      # fpga_fd
    ctypes.c_char_p,   # buffer
    ctypes.c_uint64,   # size
    ctypes.c_uint64,   # addr
]

libdma.read_to_buffer.restype = ctypes.c_ssize_t
libdma.write_from_buffer.restype = ctypes.c_ssize_t
    

# FPGA write
def richelot_isogeny_fpga(device_w, device_r, fpga_w_fd, fpga_r_fd, ker):
    write_data = b"richelot"
    write_buf = ctypes.create_string_buffer(write_data)
    rc_write = libdma.write_from_buffer(
        device_w.encode('utf-8'),
        fpga_w_fd,
        write_buf,
        len(write_data),
        0,
    )
    if rc_write < 0:
        raise RuntimeError("FPGA write failed")

       # FPGA read
    read_buf = ctypes.create_string_buffer(4096)
    rc_read = libdma.read_to_buffer(
        device_r.encode('utf-8'),
        fpga_r_fd,
        read_buf,
        128,
        ctypes.c_uint64(0),
    )
    if rc_read < 0:
        raise RuntimeError("FPGA read failed")
    print("rc_read", rc_read)
    fio.print_buffer(read_buf, rc_read)

def print_input_data(ker, h):
    D1, D2 = ker
    fio.print_polynomial_elements(D1[0], "D1[0] (f1)")
    fio.print_polynomial_elements(D1[1], "D1[1] (g1)")
    fio.print_polynomial_elements(D2[0], "D2[0] (f1)")
    fio.print_polynomial_elements(D2[1], "D2[1] (f2)")
    fio.print_polynomial_elements(h, "h")

def decompose_and_reconstruct(ker, h, K):
    """
    ker: (D1, D2) 2つの多項式係数
    h: 多項式
    K: 拡大体
    """
    R = PolynomialRing(K, name="x", implementation="generic")
    D1, D2 = ker

    # f: 多項式（次数2）、g: 定数項（多項式でない）
    f1, g1 = D1
    f2, g2 = D2
    
    # f1, g1, f2, g2
    f1_new = fio.decompose_and_reconstruct_polynomial(f1, K, R)
    g1_new = fio.decompose_and_reconstruct_polynomial(g1, K, R)
    f2_new = fio.decompose_and_reconstruct_polynomial(f2, K, R)
    g2_new = fio.decompose_and_reconstruct_polynomial(g2, K, R)
    ker_new = ((f1_new, g1_new), (f2_new, g2_new))
    # h
    h_new = fio.decompose_and_reconstruct_polynomial(h, K, R)

    return ker_new, h_new
