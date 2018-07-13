require "spec"

class BitTest < Test::Unit::TestCase
  dtype = Lattice::Bit

  test dtype do
    assert { dtype < Lattice::NArray }
  end

  procs = [
    [proc{ |tp, a| tp[*a] }, ""],
    [proc{ |tp, a| tp[*a][true] }, "[true]"],
    [proc{ |tp, a| tp[*a][0..-1] }, "[0..-1]"]
  ]
  procs.each do |init, ref|

    test ""
