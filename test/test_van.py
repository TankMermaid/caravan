#!/usr/bin/env python

'''
unit tests for van.py
'''

from caravan import van
import pytest

class TestParseArgs:
    def test_check_intersect(self):
        func, opts = van.parse_args(['check_intersect', 'for.fq', 'rev.fq'])
        assert opts == {'forward': 'for.fq', 'reverse': 'rev.fq'}