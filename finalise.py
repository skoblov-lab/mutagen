#! /usr/bin/env python

"""



"""

import os
import re
import json
from typing import TypeVar, Callable, Mapping, Union

import click
from fn import F

from mutagen import parser

ASSOC_PAT = re.compile('^[A-Z0-9]*(?::\d*){1,2}$')
PHANTOM_ASSOC_PAT = re.compile('^[A-Z0-9]*(?::[A-Z\-\d]*\?)$')


A = TypeVar('A')
B = TypeVar('B')
C = TypeVar('C')


def fnor(f1: Callable[[A], B], f2: Callable[[A], C], value: A) -> Union[B, C]:
    """
    Calling `fnor(f1, f2, value)` is equivalent to `f1(value) or f2(value)`
    :param f1: the first function to try
    :param f2: the second function to try
    :param value: a value both f1 and f2 can accept
    :return:
    """
    return f1(value) or f2(value)


def rename_associations(mapping: Mapping[str, str], record: parser.Record) \
        -> parser.Record:
    """
    Rename associations in all subrecords.
    :param mapping:
    :param record:
    :return:
    :raises ValueError: if there are missing/malformed entries at any level,
    if the mapping is malformed or non-exhaustive.
    """
    is_association = F(fnor, ASSOC_PAT.match, PHANTOM_ASSOC_PAT.match)

    if not (F(map, is_association) >> all)(mapping.values()):
        raise ValueError(f'malformed association mapping for {record.protein}')

    def check_mut(mut: parser.Mutation) -> parser.Mutation:
        if not (mut.subrecs and all(mut.subrecs)):
            raise ValueError(
                f'malformed subrecord(s) in {record.protein}:{mut.id}'
            )
        return parser.Mutation(
            mut.id, mut.start, mut.stop, mut.ref, mut.alt, mut.description,
            list(map(F(check_subrec, mut.id), mut.subrecs))
        )

    def check_subrec(mutid: int, subrec: parser.SubRecord) -> parser.SubRecord:
        if not (subrec.effects and all(subrec.effects)):
            raise ValueError(
                f'malformed effect(s) in {record.protein}:{mutid}:{subrec.id}'
            )
        effects = list(map(rename_assoc, subrec.effects))
        return parser.SubRecord(subrec.id, subrec.description, effects)

    def rename_assoc(effect: parser.Effect) \
            -> parser.Effect:
        try:
            associations = [mapping[assoc] for assoc in effect.associations]
            return parser.Effect(effect.cls, effect.level, effect.target,
                                 associations)
        except KeyError as err:
            raise ValueError(
                f'non-exhaustive association mapping in {record.protein}; {err}'
            )

    if not (record.mutations and all(record.mutations)):
        raise ValueError(f'malformed mutations in {record.protein}')

    return parser.Record(record.protein, list(map(check_mut, record.mutations)))


def optval(validator: Callable[[A], bool], message: str, ctx, param: str,
           value: A):
    if not validator(value):
        raise click.BadParameter(message, ctx=ctx, param=param)
    return value


@click.command('finalise', help=__doc__,
               context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', required=True,
              type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.option('-a', '--associations',
              type=click.Path(exists=True, dir_okay=False, resolve_path=True))
@click.option('-o', '--output', required=True,
              type=click.Path(exists=False, dir_okay=False, resolve_path=True),
              callback=F(optval, lambda v: not os.path.exists(v), 'output exists'))
def finalise(input, associations, output):
    try:
        with open(associations) as buffer:
            mappings = json.load(buffer)
    except json.JSONDecodeError:
        raise click.BadParameter("can't parse association mapping")
    try:
        with open(input) as buffer:
            records = parser.parse(parser.cleanup(buffer))
    except (ValueError, AttributeError):
        raise click.BadParameter(
            "can't parse records (this usually signifies a record-level error)"
        )
    finalised_records = [
        rename_associations(mappings.get(rec.protein, {}), rec)
        for rec in records
    ]
    with open(output, 'w') as dest:
        for line in parser.write(finalised_records):
            print(line, file=dest)


if __name__ == '__main__':
    finalise()


