using NativeFileDialog

maxSimulations = 972

folder = pick_folder(pwd()*"\\results")

if folder == ""
    return
end

folder = splitdir(folder)[end]

numberOfNumbers = 0;
while true
    try
        parse(Int,folder[end-numberOfNumbers:end])
        global numberOfNumbers += 1
    catch
        break
    end


end

if numberOfNumbers == 0
    println("No number at the end of the directory name.")
    return
end

namePart = folder[1:end-numberOfNumbers];

global A = zeros(Bool,maxSimulations)

for i = 1:maxSimulations

    global A[i] = isdir("./results/"*namePart*lpad(i,2,'0'))

end

inds = findall(.!A)