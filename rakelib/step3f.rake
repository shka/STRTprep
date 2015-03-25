##
## Step 3f - Reads and the normalized reads in the qualified samples
##


file 'out/cg/reads.txt.gz' =>
     ['out/cg/samples.txt', 'out/cg/reads_all.txt.gz'] do |t|
  qualified = Hash.new
  infp = open("| grep -v ^LIBRARY #{t.source}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    qualified["#{cols[0]}.#{cols[1]}"] = ''
  end
  infp.close

  header = [true]
  outfp = open("| gzip -c > #{t.name}", 'w')
  infp = open("| gunzip -c #{t.sources[1]}")
  while line = infp.gets
    cols = line.rstrip.split(/\t/)
    if header.length == 1
      tmp = ['Gene']
      cols[1..-1].each do |colname|
        libwellid, name = colname.split(/\|/)
        header.push(qualified.key?(libwellid))
        tmp.push(colname) if qualified.key?(libwellid)
      end
      outfp.puts tmp.join("\t")
    else
      tmp = [cols[0]]
      sum = 0
      1.upto(header.length-1) do |i|
        if header[i]
          tmp.push(cols[i])
          sum += cols[i].to_i
        end
      end
      outfp.puts tmp.join("\t") if sum > 0
    end
  end
  infp.close
  outfp.close
end

file 'out/cg/nreads.txt.gz' => 'out/cg/reads.txt.gz' do |t|
  sh "R --vanilla --quiet < bin/_step3f_normalization.R > #{t.name}.log 2>&1"
end
