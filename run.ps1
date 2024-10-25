$target = -5.6104
$nfesLmt = 1000000
$runtimeLmt = 120
$Np = 300
$sequence = "ABBBBBBABBBAB"
$seed_from = 1
$seed_to = 50
$runMode = "Release"
$outputFile = "C:\Users\Mic\PycharmProjects\protein_folding\results_${Np}_p_${nfesLmt}_nfes_from_${seed_from}_to_${seed_to}.txt"

# Clear the output file if it exists
if (Test-Path $outputFile) {
    Remove-Item $outputFile
}

"S|Time|Op/s|Energy" | Out-File -FilePath $outputFile -Append
# Loop from 1 to 10, changing the seed on each iteration
for ($seed = $seed_from; $seed -le $seed_to; $seed++) {
    # Construct the command string
    #$command = "C:\Users\Mic\source\repos\protein_folding\x64\$runMode\protein_folding.exe $sequence -seed $seed -target $target -nfesLmt $nfesLmt -runtimeLmt $runtimeLmt -Np $Np"
    $command = "C:\Users\Mic\source\repos\protein_folding\protein_folding\protein_folding.exe $sequence -seed $seed -target $target -nfesLmt $nfesLmt -runtimeLmt $runtimeLmt -Np $Np"

    # Run the command and capture the output in a variable
    $output = Invoke-Expression $command
    
    # Append seed info and output to the results file
    $output | Out-File -FilePath $outputFile -Append
    
    # Optional: Display progress in the console
    Write-Host "Iteration $seed completed."
}

# Output file will contain results from all iterations
Write-Host "All results written to $outputFile."
 
