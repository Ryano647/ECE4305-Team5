see =[];
for test = 1:20
    see =[see, data(test)== demodulatedData(test)]
end
