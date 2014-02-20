#!/usr/bin/env ruby

require 'bio-faster'

abort "Specify tiles to extract with regular expression." if ARGV.length != 1

p = Regexp.new(ARGV[0])

Bio::Faster.new(:stdin).each_record(:quality => :raw) { |id, seq, qual|
  cols = id.split(/\:/)
  if p.match(cols[4])
    puts ["@#{id}", seq, '+', qual].join("\n")
  end
}
