#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'arginine'

par = Arginine::parse do
  desc "convert allrank or fixrank result to qiime format"
  opt :threshold, desc: "truncate lineage if support falls below this fraction", default: 0.50, cast: :to_f
end

LEVELS = %w{k__ p__ c__ o__ f__ g__ s__}

in_header = true
ARGF.each do |line|
  line.chomp!
  if in_header
    in_header = false if line == ""
  else
    fields = line.split(/;/)
    name = fields.shift # the first field is the name
    fields = fields.drop_while { |f| f != "+" } # throw away stuff before "+"

    # throw away the "root" fields if they are present
    if fields[1] == "Root"
      fields.shift(3)
    else
      fields.shift(1)
    end

    LEVELS.zip(fields.each_slice(2)).map do |code, slice|
      phyl, support = slice
      if phyl.nil? or support.to_f < par[:threshold]
        print [code, "; "].join("")
      else
        print [code, phyl.gsub(/"/, ''), "; "].join("")
      end
    end
    print "\n"
  end
end