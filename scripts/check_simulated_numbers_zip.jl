using NativeFileDialog

maxSimulations = 64

file = pick_file(pwd()*".\\results";filterlist="zip")

if file == ""
    return
end

file = splitdir(file)[end][1:end-4]

numberOfNumbers = 0;
while true
    try
        parse(Int,file[end-numberOfNumbers:end])
        global numberOfNumbers += 1
    catch
        break
    end


end

if numberOfNumbers == 0
    println("No number at the end of the directory name.")
    return
end

namePart = file[1:end-numberOfNumbers];

global A = zeros(Bool,maxSimulations)

for i = 1:maxSimulations

    global A[i] = isfile("./results/"*namePart*lpad(i,numberOfNumbers,'0')*".zip")

end

inds = findall(.!A)