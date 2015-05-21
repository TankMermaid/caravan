#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'arginine'
require 'bio'

class HashSplitter
  def initialize(n_buckets, output_base)
    @n_buckets = n_buckets
    @output_base = output_base

    @n_digits = (@n_buckets - 1).to_s.length
    @output_fhs = Hash.new
  end

  def hash(entry)
    entry.data.gsub(/\n/, '').hash % @n_buckets
  end

  def output_fn(bucket)
    @output_base + "." + "%0#{@n_digits}d" % bucket
  end

  def output_fh(entry)
    h = hash(entry)
    @output_fhs[h] = open(output_fn(h), 'w') unless @output_fhs.include? h
    @output_fhs[h]
  end

  def dump(entry)
    output_fh(entry).write(entry)
  end
end

par = Arginine::parse do
  desc "split identical entries into the same buckets"
  arg :n, cast: :to_i, desc: "number of buckets"
  arg :output, desc: "output base"
  argf "fasta"
end

hs = HashSplitter.new(par[:n], par[:output])

Bio::FlatFile.foreach(ARGF) do |entry|
  hs.dump(entry)
end