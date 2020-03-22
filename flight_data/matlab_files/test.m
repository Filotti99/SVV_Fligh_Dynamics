flightdata = load("FTISxprt-20200310_flight2.mat");
flightdata_array = struct2array(flightdata.flightdata);
for flightdata_item = flightdata_array
    name = flightdata_item.description+"["+flightdata_item.units+"].csv";
    new_name = regexprep(name, "/", "_p_");
    csvwrite(new_name, flightdata_item.data);
end