#!/usr/bin/env ruby
## plot depolarization vs (x,Q2)
## run after analysis_depolarization.C -> postprocess_depolarization.C

require 'pycall/import'
r = PyCall.import_module 'ROOT'
r.gROOT.SetBatch true
r.gStyle.SetOptStat 0

if ARGV.length<1
  puts "USAGE: #{$0} [out/___.canvas.root file]"
  exit 2
end
rootFileName = ARGV[0]

# read plots
rootFile = r.TFile.new rootFileName
plotList = PyCall.iterable(rootFile.GetListOfKeys)
  .select{ |key| key.GetName.match? /^histos/ }
  .select{ |key| key.GetName.match? /Q2vsX/ }
  .map &:ReadObj

# canvas
canv = r.TCanvas.new 'depol', 'depol', 2500, 1800
canv.Divide 3, 3
plotList.each_with_index do |plot,i|

  # get maximum
  max = (1..plot.GetNbinsX).map do |bx|
    (1..plot.GetNbinsY).map do |by|
      plot.GetBinContent bx, by
    end
  end.flatten.max
  plot.GetZaxis.SetRangeUser 0, 1.5*max

  # draw
  pad = canv.GetPad i+1
  pad.cd
  pad.SetLogx
  pad.SetLogy
  pad.SetPhi   40
  pad.SetTheta 15
  pad.SetLeftMargin 0.25
  [ plot.GetXaxis, plot.GetYaxis, plot.GetZaxis ].each do |ax|
    ax.SetTitleOffset 1.3
    ax.SetTitleSize   0.05
  end
  plot.SetTitle ''
  plot.Draw 'lego2'
end

canv.SaveAs rootFileName.gsub(/root$/,"Q2vsX.png")
