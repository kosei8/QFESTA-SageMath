import ctypes
from sage.all import Polynomial

# input data functions       
def print_polynomial_elements(coeffs, label):
    coeffs_list = [decompose_element(coeffs[i]) for i in range(coeffs.degree() + 1)]
    print_coeffs(coeffs_list, label)

def print_coeffs(coeffs, label):
    print(f"{label}")
    for i, (imag, real) in enumerate(coeffs):
        print(f"deg {i}: imag: {imag}*i")
        print(f"deg {i}: real: {real}")

def decompose_and_reconstruct_polynomial(coeffs, K, R):
    """
    coeffs: (f, g) 2つの多項式係数
    K: 拡大体
    R: 多項式リング
    """
    coeffs_list = [decompose_element(coeffs[i]) for i in range(coeffs.degree() + 1)]
    coeffs_new = reconstruct_poly_from_coeffs(coeffs_list, K, R)
    return coeffs_new

def decompose_element(c):
    if isinstance(c, Polynomial):
        # 多項式ではないが多項式として扱う
        return (c[0], c[1])
    elif hasattr(c, "polynomial"):
        # 拡大体の元 → i に関する係数で展開
        poly_c = c.polynomial()
        real = poly_c[0]
        imag = poly_c[1] if poly_c.degree() >= 1 else 0
        return (imag, real)
    else:
        # 定数なら虚部なし
        return (0, c)
    
def reconstruct_poly_from_coeffs(coeffs, K, R):
    """
    coeffs: [(imag0, real0), (imag1, real1), ...]
    拡大体 K 上の多項式リング R.<x> で多項式を再構成する。
    
    たとえば coefficient d のとき:
      "a + b*i" = K(real) + K(imag)*K.gen()
    で x^d を掛け合わせる。
    """
    x = R.gen()
    poly = R(0)
    for d, (imag, real) in enumerate(coeffs):
        coef = K(real) + K(imag)*K.gen()  # a + b*i
        poly += coef * x**d
    return poly

# FPGA driver debug functions
def print_buffer(buffer: ctypes.Array, data_size: int) -> None:
    for line_start in range(0, data_size, 16):  # 16バイトずつ表示
        hex_part = ""
        ascii_part = ""

        for i in range(line_start, min(line_start + 16, data_size)):
            val = buffer[i]
            if isinstance(val, int):
                val_int = val
            else:
                val_int = ord(val)

            hex_part += f"{val_int:02x} "
            ascii_part += chr(val_int) if 32 <= val_int <= 126 else '.'

        print(f"{line_start:08x}  {hex_part:<48}  {ascii_part}")
