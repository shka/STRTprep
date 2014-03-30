#!/usr/bin/env ruby

require 'parallel'

str = ARGV[0]
orgbam = ARGV[1]
strbam = ARGV[2]

acc2dups = Hash.new
open("| samtools view #{orgbam} | cut -f 1 | sort -S #{sprintf('%d', 100/Parallel.processor_count+2)}% | uniq -c").each { |line|
  dups, acc = line.strip.split(/\s/); acc2dups[acc] = '' if dups != '1'
}

open("| bamToBed -i #{strbam}").each { |line|
  cols = line.rstrip.split(/\t/)
  puts [cols[0],
        str == 'fwd' ? cols[1] : cols[2].to_i-1,
        str == 'fwd' ? cols[1].to_i+1 : cols[2], cols[3], cols[4] ].join("\t") unless acc2dups.key?(cols[3])
}
