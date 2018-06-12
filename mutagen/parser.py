from typing import Iterator, Iterable, NamedTuple, List, Callable, Optional
from itertools import chain, groupby, filterfalse, repeat
import operator as op
import re

from fn import F


CLASSES = frozenset(['ORG', 'CEL', 'PAT', 'PRO', 'INT', 'IND', 'ENZ', 'MIM',
                     'LOC', 'TRA', 'CHA', 'CAR', 'UNK'])
LEVELS = frozenset(['--', '-', '0', '+', '++', 'r'])

MUTPREF = re.compile('^<(\d+)>')
SUBPREF = re.compile('^\[(\d+)\]')
EFFREF = re.compile('^>>')
PREF = frozenset(['<', '[', '>'])
MUTPATT = re.compile(
    '<(?P<id>\w+)>\s*(?P<ref>\w+)\|(?P<start>\d+)-(?P<stop>\d+)\|(?P<sub>\w+)\s+(?P<text>.+)'
)
SUBRECPATT = re.compile('\[(?P<id>\d+)\]\s*(?P<text>.+)')


Effect = NamedTuple('Effect', [
    ('cls', str),
    ('level', Optional[str]),
    ('target', Optional[str]),
    ('associations', List[str])
])

SubRecord = NamedTuple('SubRecord', [
    ('id', int),
    ('description', str),
    ('effects', List[Optional[Effect]])
])

Mutation = NamedTuple('Mutation', [
    ('id', int),
    ('start', int),
    ('stop', int),
    ('ref', Optional[str]),
    ('alt', Optional[str]),
    ('description', str),
    ('subrecs', List[Optional[SubRecord]])
])

Record = NamedTuple('Record', [
    ('protein', str),
    ('mutations', List[Optional[Mutation]])
])


def optionable(*exceptions):
    def wrapper(f):
        def opt(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                return None
        return opt
    return wrapper


def breakby(condition, lines: Iterable[str]) -> Iterator[Iterator[str]]:
    return (
        F(filterfalse, op.itemgetter(0)) >> (map, op.itemgetter(1))
    )(groupby(lines, condition))


def separate(condition: Callable[[str], bool], sep: str,
             lines: Iterable[str]) -> Iterator[str]:
    """
    Add a separator on condition
    :param lines:
    :return:
    """
    return chain.from_iterable(
        (sep, line) if condition(line) else (line,) for line in lines
    )


def cleanup(handle: Iterable[str]) -> Iterable[str]:
    """
    Remove comments from a handle to be parsed
    """
    def start(line: str) -> Optional[int]:
        candidates = list(
        filterfalse(F(op.eq, -1), [line.find('!!!'), line.find('***')])
        )
        return min(candidates) if candidates else None

    return (line[:start(line)] for line in handle)


def parse(handle: Iterable[str]) -> List[Record]:
    def parse_rec(lines: Iterator[str]) -> Record:
        protein = next(lines)
        mutations = (
            F(separate, MUTPREF.match, '') >> (breakby, op.not_) >>
            (map, parse_mut) >> list
        )(lines)
        return Record(protein, mutations)

    @optionable(ValueError, AttributeError)
    def parse_mut(lines: Iterator[str]) -> Mutation:
        id_, ref, start, stop, alt, text = MUTPATT.match(next(lines)).groups()
        ref, alt = map(lambda x: None if x.lower() == 'none' else x, [ref, alt])
        start, stop = map(int, [start, stop])
        subrecords = (
            F(separate, SUBPREF.match, '') >> (breakby, op.not_) >>
            (map, parse_subrec) >> list
        )(lines)
        return Mutation(int(id_), start, stop, ref, alt, text.strip(), subrecords)

    @optionable(ValueError)
    def parse_subrec(lines: Iterator[str]) -> SubRecord:
        id_, text = SUBRECPATT.match(next(lines)).groups()
        effects = map(parse_effect, lines)
        return SubRecord(int(id_), text, list(effects))

    @optionable(ValueError)
    def parse_effect(line: str) -> Effect:
        # there can be at most 3 blank positions
        # TODO add validators
        cls, level, target, associations, *_ = chain(
            line.lstrip('>> ').split('|'), repeat('', 3)
        )
        if cls.upper() not in CLASSES:
            raise ValueError
        level, target = map(lambda x: x.strip('?') or None, [level, target])
        return Effect(cls, level, target,
                      list(filter(bool, associations.split(';'))))

    return (
        F(map, str.strip) >> (filter, bool) >>
        (separate, lambda l: l[0] not in PREF, '') >> (breakby, op.not_) >>
        (map, parse_rec) >> list
    )(handle)


def write(records: Iterable[Record]) -> Iterable[str]:

    def indent(level: int, text: str) -> str:
        return '\t' * level + text

    def write_record(rec: Record) -> Iterable[str]:
        yield rec.protein
        yield from (F(map, write_mut) >> chain.from_iterable)(rec.mutations)

    def write_mut(mut: Mutation) -> Iterable[str]:
        position = f'{mut.ref}|{mut.start}-{mut.stop}|{mut.alt}'
        yield indent(1, f'<{mut.id}> {position} {mut.description}')
        yield from (F(map, write_subrec) >> chain.from_iterable)(mut.subrecs)

    def write_subrec(subrec: SubRecord) -> Iterable[str]:
        yield indent(2, f'[{subrec.id}] {subrec.description}')
        yield from map(write_eff, subrec.effects)

    def write_eff(eff: Effect) -> str:
        assoc = ';'.join(eff.associations)
        effect = '|'.join([val or '' for val in eff[:3]] + [assoc])
        return indent(3, f">> {effect}")

    return (F(map, write_record) >> chain.from_iterable)(records)


if __name__ == '__main__':
    raise RuntimeError




