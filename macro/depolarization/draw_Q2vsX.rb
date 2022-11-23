#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks

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
  .reject{ |key| key.GetName.match? /epsilon/ }
  .map &:ReadObj

# canvas
nrows = 3
ncols = 3
canv = r.TCanvas.new 'depol', 'depol', nrows*1000, ncols*800
canv.SetTopMargin  0.25
canv.SetLeftMargin 0.25
canv.Divide ncols, nrows
plotList.each_with_index do |plot,i|

  # get maximum
  max = (1..plot.GetNbinsX).map do |bx|
    (1..plot.GetNbinsY).map do |by|
      plot.GetBinContent bx, by
    end
  end.flatten.max
  plot.GetZaxis.SetRangeUser 0, 1.0*max

  # get depol variable name; associate to pad number
  varName = plot.GetName.split('_').find{|tok|tok.match?(/^depol|^epsilon/)}.sub(/vs.*/,'')
  padHash = {
    :depolA  => 1,
    :depolB  => 2,
    :depolV  => 3,
    :depolBA => 5,
    :depolVA => 6,
    :depolCA => 8,
    :depolWA => 9,
    :epsilon => 12,
  }

  # function to format titles
  def formatTitle(s)
    s
      .gsub( /(epsilon)/,   'var\1(x,Q^{2},y)'  )
      .gsub( /(A|B|C|V|W)/, '\1(#varepsilon,y)' )
      .gsub( /\//,          ' / '               )
  end

  # draw
  r.gStyle.SetPalette r.kBlueRedYellow
  pad = canv.GetPad padHash[varName.to_sym].to_i
  pad.cd
  pad.SetGrid 1, 1
  pad.SetLogx
  pad.SetLogy
  pad.SetPhi   40
  pad.SetTheta 15
  pad.SetLeftMargin   0.25
  pad.SetBottomMargin 0.15
  pad.SetTopMargin    0.10
  pad.SetRightMargin  0.12
  [ plot.GetXaxis, plot.GetYaxis, plot.GetZaxis ].each do |ax|
    ax.SetTitleOffset 1.4
    ax.SetTitleSize   0.05
    ax.SetLabelSize   0.05
  end
  [ plot.GetXaxis, plot.GetYaxis ].each do |ax| ax.SetNdivisions 30 end
  plot.SetTitle ' '*8 + formatTitle(plot.GetZaxis.GetTitle)
  plot.GetZaxis.SetTitle ''
  plot.Draw 'colz'
end

# decorations
canv.cd
lines = [
  r.TLine.new(0.70, 0, 0.70, 1),
  r.TLine.new(0, 0.33, 1, 0.33),
  r.TLine.new(0, 0.67, 1, 0.67),
].each do |l|
  l.SetLineWidth 2
  l.SetLineColor r.kGray+1
  l.SetNDC
  l.Draw
end
placePol = Proc.new do |padNum,pol|
  canv.cd padNum
  t = r.TLatex.new
  t.SetTextSize 0.08
  t.DrawLatexNDC 0.0, 0.5, pol
end
placePol.call 1, 'UU'
placePol.call 4, 'UL, UT'
placePol.call 7, 'LU, LL, LT'

canv.SaveAs rootFileName.gsub(/root$/,"Q2vsX.pdf")
