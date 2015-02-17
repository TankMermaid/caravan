#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'arginine'

params = Arginine::parse do
  desc ""
  arg :level
  arg :lineages
  arg :otu_table
  opt :count_cast, default: "i", desc: "cast to integer (i) for counts or float (f) for r.a."
end

raise RuntimeError, "count_cast must be i or f" unless %w(i f).include? params[:count_cast]
count_cast = ("to_" + params[:count_cast]).to_sym

level_i = %w{k p c o f g s}.index(params[:level])
raise RuntimeError, "level #{params[:level]} not recognized" if level_i.nil?

lin_fh = open(params[:lineages])
otu_fh = open(params[:otu_table])
header = otu_fh.gets
n_samples = header.split.length - 1 # skip the otu table header

tally = Hash.new { |h, k| h[k] = [0] * n_samples }
lin_fh.each_line.zip(otu_fh.each_line) do |lin_line, otu_line|
  # read in the lineage and truncate it
  lin = lin_line.split(/;\s*/)[0..level_i].join(";")

  # read in the counts from the otu_table and add them to the right lineage
  otu_counts = otu_line.split.drop(1).map(&count_cast)
  otu_counts.each_with_index { |c, i| tally[lin][i] += c }
end

# sum up the counts for each lineage
sums = tally.each_key.inject(Hash.new) { |h, lin| h[lin] = tally[lin].reduce(:+); h }

# make a list of lineages sorted by decreasing abundance
lins = tally.keys
lins.sort_by! { |k| -sums[k] }

# output
puts header
lins.each { |lin| puts ([lin] + tally[lin]).join("\t") }