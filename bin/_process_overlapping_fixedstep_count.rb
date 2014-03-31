#!/usr/bin/env ruby

require 'parallel'

def report_count(chr, pos, str, cnt)
  pchr = chr =~ /RNA_SPIKE_/ ? (pos == 1 && str == '+' ? chr : chr.downcase) : chr
  puts "#{pchr}:#{pos}-#{pos+49},#{str}\t#{cnt}"
end

[0, 1].each { |i|
  chr = ''
  str = i == 0 ? '+' : '-'
  pos1 = pos2 = cnt1 = cnt2 = 0
  open("| gunzip -c #{ARGV[i]} | cut -f 1,3 | sort -k 1,1 -k 2,2n -S #{sprintf('%d', 100/(Parallel.processor_count+2))}% | uniq -c").each { |line|
    cnts, tmpchr, poss = line.strip.split(/\s/)
    if tmpchr != chr
      if chr != ''
        report_count(chr, pos1, str, cnt1) if 0 < pos1
        report_count(chr, pos2, str, cnt2) if 0 < pos2
      end
      chr = tmpchr
      pos1 = pos2 = cnt1 = cnt2 = 0
    end
    pos = poss.to_i
    cnt = cnts.to_i
    if pos1 == 0
      pos1 = (pos-1)/50*50+1
      cnt1 = cnt
    elsif pos1+49 < pos
      report_count(chr, pos1, str, cnt1)
      pos1 = (pos-1)/50*50+1
      cnt1 = cnt
    else
      cnt1 += cnt
    end
    if pos2 == 0
      if 25 < pos
        pos2 = (pos-26)/50*50+26
        cnt2 = cnt
      end
    elsif pos2+49 < pos
      report_count(chr, pos2, str, cnt2)
      pos2 = (pos-26)/50*50+26
      cnt2 = cnt
    else
      cnt2 += cnt
    end
  }
  if chr != ''
    report_count(chr, pos1, str, cnt1) if 0 < pos1
    report_count(chr, pos2, str, cnt2) if 0 < pos2
  end
}
