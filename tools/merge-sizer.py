#!/usr/bin/env python3

import textwrap, curses, argparse
from Bio import SeqIO

class Alignment:
    def __init__(self, win, forward, reverse, offset=0):
        self.win = win
        self.offset = offset

        self.forward = forward
        self.reverse = reverse
        self.forward_records = SeqIO.parse(forward, 'fastq')
        self.reverse_records = SeqIO.parse(reverse, 'fastq')

        self.max_y, self.max_x = win.getmaxyx()
        self.margin = 2
        self.wrap_width = self.max_x - self.margin

        self.next_record()

    def next_record(self):
        self.forward_record = next(self.forward_records)
        self.reverse_record = next(self.reverse_records)
        self.forward_seq = str(self.forward_record.seq)
        self.reverse_seq = str(self.reverse_record.reverse_complement().seq)
        self.display()
        
    @staticmethod
    def offset_strings(s1, s2, offset):
        empties = "-" * abs(offset)
        if offset == 0:
            return (s1, s2)
        elif offset > 0:
            return (s1 + empties, empties + s2)
        else:
            return (empties + s1, s2 + empties)

    @staticmethod
    def alignment_char(c1, c2):
        if c1 == c2:
            return "|"
        else:
            return " "

    @classmethod
    def alignment(cls, s1, s2):
        return "".join([cls.alignment_char(c1, c2) for c1, c2 in zip(s1, s2)])

    def change_offset_by(self, change): 
        self.offset += change
        self.display()

    def change_offset_to(self, offset):
        self.offset = offset
        self.display()

    @classmethod
    def wrap(cls, s, width):
        if len(s) < width:
            return [s]
        else:
            return [s[0:width]] + cls.wrap(s[width:], width)


    def display(self):
        # compute offset strings
        t1, t2 = self.offset_strings(self.forward_seq, self.reverse_seq, self.offset)
        aln = self.alignment(t1, t2)
        lines_groups = zip(*[self.wrap(x, width=self.wrap_width) for x in [t1, aln, t2]])

        self.win.erase()
        y = 1
        for lines_group in lines_groups:
            for line in lines_group:
                self.win.addstr(y, 0, line)
                y += 1
            y += 1

        y += 1
        y += 1
        self.win.addstr(y, 0, self.forward_record.id)
        y += 1
        self.win.addstr(y, 0, self.reverse_record.id)
        y += 1
        self.win.addstr(y, 0, "offset: {}".format(self.offset))
        y += 1
        self.win.addstr(y, 0, "trimmed construct size: {}".format(sum([1 for a, b in zip(t1, t2) if a != "-" and b != "-"])))
        y += 1
        self.win.addstr(y, 0, "total construct size: {}".format(max(len(t1), len(t2))))
        y += 1
        self.win.addstr(y, 0, "matching positions: {}".format(aln.count("|")))
        y += 1
        self.win.addstr(y, 0, "unmatched positions: {}".format(len(aln) - aln.count("|")))

        self.win.refresh()


def keyloop(stdscr, forward, reverse):
    stdscr.clear()
    stdscr_y, stdscr_x = stdscr.getmaxyx()
    subwin = stdscr.subwin(stdscr_y - 3, stdscr_x, 0, 0)
    align = Alignment(subwin, forward, reverse)
    
    while True:
        c = stdscr.getch()
        if 0 < c < 256:
            c = chr(c)
            if c in 'Qq':
                break
            elif c in 'Hh':
                align.change_offset_by(-1)
            elif c in 'Ll':
                align.change_offset_by(1)
            elif c in 'Jj':
                align.next_record()
            elif c == '0':
                align.change_offset_to(0)
            else:
                pass
        else:
            if c == curses.KEY_LEFT:
                align.change_offset_by(-1)
            elif c == curses.KEY_RIGHT:
                align.change_offset_by(1)
            elif c == curses.KEY_DOWN:
                align.next_record()
            else:
                pass

def main(stdscr, forward, reverse):
    keyloop(stdscr, forward, reverse)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('forward', type=argparse.FileType('r'))
    p.add_argument('reverse', type=argparse.FileType('r'))
    args = p.parse_args()

    curses.wrapper(main, args.forward, args.reverse)
