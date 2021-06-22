#migration p1~p2 Scenario 2
python3 wrapper_try.py migration "" 3 "100000,50000" "" "0,1/1,0=0.0002/0.002/0.05/0.2/0.5" ""

#migration and bottleneck Scenario 4
python3 wrapper_try.py migration bottleneck 3 "100000,50000" "" "0,1/1,0=0.0002/0.002/0.05/0.2" "ppc:2000:500:0/ppc:10000-50000-10000:10000:0"
python3 wrapper_try.py migration bottleneck 3 "100000,50000" "" "0,1/1,0=0.0002/0.002/0.05/0.2" "ppc:2000:500:1/ppc:10000-50000-10000:10000:1"

#migration p1~p2 Scenario 6
python3 wrapper_try.py migration "" 3 "100000,50000" "" "0,1/1,0=0.000007/0.00007/0.00009/0.0001/0.0003/0.0005/0.0007" ""


