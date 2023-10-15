ENV["LUMERICAL_PYTHON_PATH"]="/opt/lumerical/v212/api/python/"
using NanophotonicModeling_LumericalFDTD,NanophotonicModeling
material_defs=Dict("TiN"=>"TiN - Palik","Si"=>"Si (Silicon) - Palik","air"=>"etch","SiO2"=>"SiO2 (Glass) - Palik")
l1=PatternedLayer(100,"TiN",[Inclusion("air",Circle(500.0))])
l2=PlainLayer(100,"Si")
l3=PlainLayer(2000,"SiO2")
mdl=PeriodicLayerModel([l1,l2,l3],"air","Si",[1000.0,1000])

hp::Dict{Symbol,Any}=Dict(:material_defs=>material_defs)

NanophotonicModeling_LumericalFDTD.build_model(mdl,hp,"/tmp/test.fsp")
