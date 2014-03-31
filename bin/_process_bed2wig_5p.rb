#!/usr/bin/env ruby

require 'parallel'

libwellid = ARGV[0]
str = ARGV[1]
bed = ARGV[2]

offset = str == 'fwd' ? 1 : 0
signal = str == 'fwd' ? 1 : -1
chr2pos2cnt = Hash.new
open("| gunzip -c #{bed} | cut -f 1,#{str == 'fwd' ? 2 : 3},4").each { |line|
  chr, poss, acc = line.rstrip.split(/\t/)
  next unless chr =~ /^chr/
  if !chr2pos2cnt.key?(chr)
    chr2pos2cnt[chr] = Hash.new
  end
  pos = poss.to_i + offset
  if !chr2pos2cnt[chr].key?(pos)
    chr2pos2cnt[chr][pos] = signal
  else
    chr2pos2cnt[chr][pos] += signal
  end
}

puts "track type=wiggle_0 name=\"#{libwellid}\" description=\"\" alwayzZero=on visivility=dense maxHeightPixel=64:32:16 color=#{str != 'fwd' ? '255,128,0' : '0,128,255'}"
chr2pos2cnt.each { |chr, pos2cnt|
  puts "variableStep chrom=#{chr}"
  pos2cnt.keys.sort.each { |pos|
    puts [pos, pos2cnt[pos]].join("\t")
  }
}
