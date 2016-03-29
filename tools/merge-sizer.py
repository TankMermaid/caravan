#!/usr/bin/env python3

# author: Scott W. Olesen <swo@mit.edu>
# to do: add a menu, error flash, record number, help screen

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

        self.record_i = -1
        self.record_buffer = []
        self.next_record()

    def next_record(self):
        if self.record_i + 1 >= len(self.record_buffer):
            # need to get a new record
            try:
                forward_record = next(self.forward_records)
                reverse_record = next(self.reverse_records)

                self.record_buffer.append((forward_record, reverse_record))
                self.record_i += 1
                assert self.record_buffer[self.record_i] == (forward_record, reverse_record)
            except StopIteration:
                pass
        else:
            self.record_i += 1

        self.update_record()

    def previous_record(self):
        if self.record_i > 0:
            self.record_i -= 1

        self.update_record()

    def update_record(self):
        self.forward_record, self.reverse_record = self.record_buffer[self.record_i]

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
        if c1 == c2 and c1 != "-":
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
        self.win.addstr(y, 0, "record #{}: id={}".format(self.record_i + 1, self.forward_record.id))
        y += 1
        self.win.addstr(y, 0, "        {}  id={}".format(" " * len(str(self.record_i + 1)), self.reverse_record.id))
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


class Pender:
    def __init__(self, align):
        self.mode = 'command'
        self.pending = []
        self.align = align

        self.operator_keys = ['g', 'h', 'j', 'k', 'l', curses.KEY_LEFT, curses.KEY_RIGHT, curses.KEY_UP, curses.KEY_DOWN]

    def send_char(self, c):
        if 0 < c < 256:
            c = chr(c)

        if self.mode == 'command' and c in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '-']:
            self.mode = 'pending'
            self.pending = [c]
        elif self.mode == 'pending' and c in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            self.pending.append(c)
        elif self.mode == 'pending' and c == 27: 
            # escape
            self.mode = 'command'
        elif self.mode == 'pending' and c in self.operator_keys:
            magnitude = int("".join(self.pending))
            if c == 'g':
                self.align.change_offset_to(magnitude)
            elif c in ['h', curses.KEY_LEFT]:
                self.align.change_offset_by(-magnitude)
            elif c in ['l', curses.KEY_RIGHT]:
                self.align.change_offset_by(magnitude)
            elif c in ['j', curses.KEY_DOWN]:
                for i in range(magnitude):
                    self.align.next_record()
            elif c in ['k', curses.KEY_UP]:
                for i in range(magnitude):
                    self.align.previous_record()

            self.mode = 'command'
        elif self.mode == 'command':
            if c == '0':
                self.align.change_offset_to(0)
            elif c in ['j', curses.KEY_DOWN]:
                self.align.next_record()
            elif c in ['k', curses.KEY_UP]:
                self.align.previous_record()
            elif c in ['h', curses.KEY_LEFT]:
                self.align.change_offset_by(-1)
            elif c in ['l', curses.KEY_RIGHT]:
                self.align.change_offset_by(1)
            elif c == 'q':
                # break if signalled by returning false
                return False
            else:
                pass
        else:
           # should have an error saying i don't know what you pressed
           self.mode = 'command'

        return True


def keyloop(stdscr, forward, reverse):
    stdscr.clear()
    stdscr_y, stdscr_x = stdscr.getmaxyx()
    subwin = stdscr.subwin(stdscr_y - 3, stdscr_x, 0, 0)
    align = Alignment(subwin, forward, reverse)
    pend = Pender(align)
    
    while True:
        c = stdscr.getch()
        status = pend.send_char(c)

        if not status:
            break

def main(stdscr, forward, reverse):
    keyloop(stdscr, forward, reverse)

if __name__ == '__main__':
    p = argparse.ArgumentParser(
    description='''
This script should help you look at some unmerged fastq's and decide on whether
the reads are staggered and what size of merged construct you can expect.
    
The script moves through pairs of entries in the two fastq's. The forward read
will be shown on the top of the two lines; the (reverse completement of) the
reverse read on the lower line. An alignment marker (showing '|' if the
corresponding bases are equal) is in between. There are some output showing
where you are in the fastq's, the offset you chose, size of the construct, etc.
   
Commands (vim-like):
    j or up: move back an entry
    k or down: move forward an entry
    h or left: move the reverse read to the left
    l or right: move the reverse read to the right
    [-0-9] + g: go to offset specified by previous number
    [-0-9] + movement key: move left/right/up/down according to number
    esc: cancel numbers added up so far
    0: go to zero offset
    q: quit''', formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('forward', type=argparse.FileType('r'), help='forward fastq')
    p.add_argument('reverse', type=argparse.FileType('r'), help='reverse fastq')
    args = p.parse_args()

    curses.wrapper(main, args.forward, args.reverse)
