#!/usr/bin/env ruby

require 'bio-faster'

n = ARGV[0].nil? 500000 : ARGV[0].to_i

outfp = IO.try_convert(STDOUT)

pretile = 0
i = -1
out = ''
Bio::Faster.new(:stdin).each_record(:quality => :raw) do |acc, seq, qv|
  tile = acc.split(':')[4][-2..-1].to_i
  if pretile != tile
    i = n
    pretile = tile
  end
  if i > 0
    out << "@#{acc}\n#{seq}\n+\n#{qv}\n"
    while 0 < out.bytesize
      begin
        written = outfp.write_nonblock(out)
        out = out.byteslice(written..-1)
      rescue IO::WaitWritable
        if 1048576 < out.bytesize
          IO.select(nil, [outfp])
          retry
        end
      end
    end    
    i -= 1
  end
end
outfp.write(out) if 0 < out.bytesize
outfp.close
