from dataclasses import dataclass
import logging
import pyarrow
import pyarrow.csv
import regex
from phantombuster.plumbing import calculate_threshold

VALID_TAGS = ['b2', 'query', 'bc', 'name', 'seq']


def read_file(name: str) -> pyarrow.Table:
    if name.endswith('tsv'):
        options = pyarrow.csv.ParseOptions(delimiter="\t")
    else:
        options = pyarrow.csv.ParseOptions()
    table = pyarrow.csv.read_csv(name, parse_options=options)
    return table


@dataclass
class InputFile:
    path: str
    group: str
    prefix: str

    def _to_json_e(self):
        return {'path': self.path, 'group': self.group, 'prefix': self.prefix}

    @classmethod
    def _from_json_e(cls, d: dict):
        return InputFile(**d)


@dataclass
class InputGroup:
    files: list[InputFile]
    group: str

    def _to_json_e(self) -> dict:
        return {'files': [f._to_json_e() for f in self.files], 'group': self.group}

    @classmethod
    def _from_json_e(cls, d: dict) -> 'InputGroup':
        ig = InputGroup(files=[InputFile._from_json_e(f) for f in d['files']], group=d['group'])
        return ig


class RegexDictionary:

    def __init__(self, d=None):
        if d is None:
            d = {}
        self._groups = d

    def __eq__(self, o):
        if isinstance(o, RegexDictionary):
            return self._groups == o._groups
        else:
            return NotImplemented()

    def add_regex(self, tag, rex, prefix='', group='*'):
        d = self._groups.get(group, {})
        d[prefix+tag] = rex
        self._groups[group] = d

    def get_regexes_for_group(self, group):
        d1 = self._groups.get(group, {})
        d2 = self._groups.get('*', {})
        d = d1|d2
        assert len(d) > 0
        return {key: regex.compile(re) for key, re in d.items()}

    def _to_json_e(self):
        return self._groups

    @classmethod
    def _from_json_e(cls, d):
        r = RegexDictionary(d)
        return r


def read_regex_file(path):
    table = read_file(path)
    column_names = ['group', 'prefix', 'tag', 'regex']

    tags = table['tag'].to_pylist()

    invalid_tags = [tag for tag in tags if tag not in VALID_TAGS]
    if len(invalid_tags) > 0:
        raise KeyError(f'Invalid tags: {", ".join(invalid_tags)}')

    regexes = table['regex'].to_pylist()

    if 'group' not in table.column_names:
        groups = ['*'] * len(tags)
    else:
        groups = table['group'].to_pylist()

    if 'prefix' not in table.column_names:
        prefixs = [''] * len(tags)
    else:
        prefixs = table['prefix'].to_pylist()

    regex_dict = RegexDictionary()

    for group, tag, re, prefix in zip(groups, tags, regexes, prefixs):
        regex_dict.add_regex(tag, re, prefix, group)
    return regex_dict


def read_barcode_hierarchy_file(name):
    table = read_file(name)
    column_names = ["barcode", "type", "referencefile", "threshold", "min_length", "max_length"]
    missing_columns = [col for col in column_names if col not in table.column_names]
    if len(missing_columns) > 0 :
        raise Exception(f"Barcode Hierarchy File Incorrect, missing columns: {missing_columns}") # TODO better error reporting

    def parse_length(value):
        if value == "-":
            return None
        else:
            try:
                return int(value)
            except Exception:
                raise Exception("Can not parse min or max length correctly in barcode hierarchy file, must be either '-' or an integer")

    def parse_type(value):
        if value == "reference":
            return "reference"
        else:
            return "random"

    barcodes = [{"name": name, "type": parse_type(type), "referencefile": reference, "threshold": threshold,
                 "min_length": parse_length(min_length), "max_length": parse_length(max_length)}
                 for name, type, reference, threshold, min_length, max_length in zip(*[table[name].to_pylist() for name in column_names])]

    for bc in barcodes:
        if bc["type"] == "reference":
            table = read_file(bc["referencefile"])
            bc["reference"] = {bc: name for bc, name in zip([str(bc) for bc in table["barcode"].to_pylist()], [str(name) for name in table["name"].to_pylist()])}
            if bc["threshold"] == "auto":
                bc["threshold"] = calculate_threshold(bc["reference"])
            else:
                try:
                    bc["threshold"] = int(bc["threshold"])
                except Exception:
                    raise Exception("Can not properly read threshold column of barcode hierarchy file, needs to be 'auto' or an integer")
    return barcodes


def read_input_files_file(path):
    table = read_file(path)

    column_names = ['file', 'group', 'prefix']

    files = table['file'].to_pylist()

    if 'group' not in table.column_names:
        groups = table["file"].to_pylist()
    else:
        groups = table['group'].to_pylist()

    if 'prefix' not in table.column_names:
        prefixs = [''] * len(files)
    else:
        prefixs = table['prefix'].to_pylist()

    input_files = [InputFile(file, group if group is not None else "", prefix if prefix is not None else "") for file, group, prefix in zip(files, groups, prefixs)]
    groups = {}
    for f in input_files:
        l = groups.get(f.group, [])
        l.append(f)
        groups[f.group] = l

    igs = []
    for group_name, files in groups.items():
        igs.append(InputGroup(files, group_name))

    return igs
