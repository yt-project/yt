"""
Static read/write methods for idefix .ini files.
"""
import re

from yt.funcs import ensure_list

_section_exp = re.compile(r"\[\w+\]\s*")
_sci_notation_exp = re.compile(r"\d+(\.\d*)?e[+-]?\d+?")


def _decode_sci_int(s):
    """
    Cast an 'e' formatted string `s` to integer if such a conversion can
    be perfomed without loss of data. Raise ValueError otherwise.

    Examples
    --------
    >>> _decode_sci_int("6.28E2")
    628
    >>> _decode_sci_int("1.4e3")
    1400
    >>> _decode_sci_int("7.0000E2")
    700
    >>> _decode_sci_int("700.00E-2")
    7
    >>> _decode_sci_int("700e-2")
    7
    >>> _decode_sci_int("6.000e0")
    6
    >>> _decode_sci_int("700e-3")
    Traceback (most recent call last):
    ...
    ValueError
    >>> _decode_sci_int("7.0001E2")
    Traceback (most recent call last):
    ...
    ValueError
    >>> _decode_sci_int("0.6e0")
    Traceback (most recent call last):
    ...
    ValueError

    """
    s = s.lower()  # assuming Idefix knows how to read "1e3" as well as "1E3"

    if not re.match(_sci_notation_exp, s):
        raise ValueError

    digits, exponent = s.split("e")
    exponent = int(exponent)
    if "." in digits:
        digits, decimals = digits.split(".")
        decimals = decimals.rstrip("0")
    else:
        decimals = ""

    if exponent >= 0 and len(decimals) > exponent:
        raise ValueError
    elif exponent < 0 and len(digits) - 1 < -exponent:
        raise ValueError

    return int(float(s))


def _encode_sci(r):
    """
    Convert a real number `r` to string, using scientific notation.

    Note that this differs from using format specifiers (e.g. `.6e`)
    in that trailing zeros are removed.
    Precision must be conserved.

    Parameters
    ----------
    r: real number (float or int)

    Returns
    -------
    ret: str
        A string representing a number in sci notation.

    Examples
    --------
    >>> _encode_sci(1)
    '1e0'
    >>> _encode_sci(0.0000001)
    '1e-7'
    >>> _encode_sci(10_000_000)
    '1e7'
    >>> _encode_sci(156_000)
    '1.56e5'
    >>> _encode_sci(0.0056)
    '5.6e-3'
    >>> _encode_sci(3.141592653589793)
    '3.141592653589793e0'
    """
    max_ndigit = len(str(r).replace(".", "")) - 1
    fmt = f".{max_ndigit}e"
    s = "{:^{}}".format(r, fmt).replace("+", "")
    ret = re.sub(r"\.?0*(e-?)0", r"\1", s)
    if ret.endswith("-"):
        ret = ret.replace("-", "0")
    return ret


def _encode_preferential_sci(r):
    """
    Convert a real number `r` to string, using sci notation if
    and only if it saves space.

    Examples
    --------
    >>> _encode_preferential_sci(189_000_000)
    '1.89e8'
    >>> _encode_preferential_sci(189)
    '189'
    >>> _encode_preferential_sci(900)
    '900'
    >>> _encode_preferential_sci(1)
    '1'
    >>> _encode_preferential_sci(0.7)
    '0.7'
    >>> _encode_preferential_sci(0.00007)
    '7e-5'
    """
    return min(str(r), _encode_sci(r), key=lambda x: len(x))


class IdefixConf(dict):
    def __init__(self, dict_or_path_or_buffer):
        """
        Parse a .ini idefix configuration from a file, or a dict.

        Parameters
        ----------
        dict_filepath_or_buffer: dict, os.Pathlike or str or a readable file handle.
        """
        if isinstance(dict_or_path_or_buffer, dict):
            super(IdefixConf, self).__init__(dict_or_path_or_buffer)
            return
        self.from_file(dict_or_path_or_buffer)

    def from_file(self, filepath_or_buffer):
        _dict = {}
        try:
            data = filepath_or_buffer.read()
        except AttributeError:
            # this is a path
            with open(filepath_or_buffer, mode="rt") as fh:
                data = fh.read()
        lines = IdefixConf.normalize_data(data)

        for line in lines:
            section = re.match(_section_exp, line)
            if section is not None:
                current_section = section.group().strip("[]")
                _dict[current_section] = {}
                continue

            key, values = IdefixConf.tokenize_line(line)
            if len(values) == 1:
                values = values[0]
            _dict[current_section][key] = values
        super(IdefixConf, self).__init__(_dict)

    @staticmethod
    def normalize_data(data):
        # normalize text body `data` to parsable text lines
        lines = []
        for line in data.split("\n"):
            if "#" in line:
                # remove comments
                line = line[: line.index("#")]

            # normalize whitespace
            line = line.strip()
            line = re.sub(r"\s", " ", line)
            if line != "":
                lines.append(line)
        return lines

    @staticmethod
    def tokenize_line(line):
        key, *raw_values = line.split()
        if not raw_values:
            raise ValueError(f"Could not parse invalid line\n{line}")

        values = []
        for val in raw_values:
            # remove period and trailing zeros to cast to int when possible
            val = re.sub(r"\.0+$", "", val)
            for caster in [int, _decode_sci_int, float, str]:
                # cast to types from stricter to most permissive
                # `str` will always succeed since it is the input type
                try:
                    values.append(caster(val))
                    break
                except ValueError:
                    continue

        return key, values

    def write_to_buffer(self, buffer):
        is_first = True
        for section, data in self.items():
            if not is_first:
                buffer.write("\n\n")
            lines = []
            buffer.write(f"[{section}]\n\n")
            for key, val in data.items():
                line = f"{key}  "
                str_val = []
                for v in ensure_list(val):
                    if isinstance(v, (float, int)):
                        str_v = _encode_preferential_sci(v)
                    else:
                        str_v = str(v)
                    str_val.append(str_v)
                val = "  ".join([v for v in str_val])
                line += str(val)
                lines.append(line)

            buffer.write("\n".join(lines) + "\n")
            is_first = False

    def write_to_file(self, filepath):
        with open(filepath, mode="wt") as fh:
            self.write_to_buffer(fh)

    def write(self, filepath_or_buffer):
        """
        Dump a dict to file in idefix's .ini format.
        """
        try:
            self.write_to_buffer(filepath_or_buffer)
        except AttributeError:
            self.write_to_file(filepath_or_buffer)
