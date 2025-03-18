#!/bin/bash

# Input file
# Uncomment the desired input file
input_file=$(ls gisaid_auspice_input_hcov-19*/*metadata.tsv 2>/dev/null | head -n 1)

# Output files
date_file="dates.tsv"
dates_decimal_file="dates_decimal.tsv"
location_metadata="location_metadata.tsv"
location_file="location_file.tsv"

# Extract date metadata (strain and date columns)
awk -F'\t' 'NR==1 {print $1 "\t" $4} NR>1 {print $1 "\t" $4}' "$input_file" > "$date_file"

# Extract location metadata (strain, country, division, and location columns)
awk -F'\t' 'NR==1 {print $1 "\t" $6 "\t" $7 "\t" $8} NR>1 {print $1 "\t" $6 "\t" $7 "\t" $8}' "$input_file" > "$location_metadata"

# Create location_file (strain and country, with no spaces in the country column)
awk -F'\t' '
BEGIN { OFS="\t"; print "strain", "country" }
NR > 1 {
    gsub(/ /, "", $6); # Remove spaces from the country column
    print $1, $6;
}
' "$input_file" > "$location_file"

# Generate dates_decimal file
awk -F'\t' '
BEGIN {
    OFS="\t";
    print "strain", "decimal_date"
}
NR > 1 {
    split($2, d, "-"); # Adjusted to correctly split date from the input file
    year = d[1];
    month = d[2];
    day = d[3];
    if (year && month && day) { # Check if the date is valid
        # Calculate day of the year
        doy = day_of_year(year, month, day);
        # Get decimal date
        decimal_date = year + (doy / days_in_year(year));
        printf "%s\t%.8f\n", $1, decimal_date; # High precision decimal date
    } else {
        print $1, "NA"; # Handle cases with missing or invalid dates
    }
}
function days_in_year(year) {
    # Check for leap year
    return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)) ? 366 : 365;
}
function day_of_year(year, month, day) {
    # Array of days in each month
    m_days = "31 28 31 30 31 30 31 31 30 31 30 31";
    split(m_days, days, " ");
    # Adjust for leap year
    if (days_in_year(year) == 366) {
        days[2] = 29;
    }
    sum = 0;
    for (i = 1; i < month; i++) {
        sum += days[i];
    }
    return sum + day;
}
' "$date_file" > "$dates_decimal_file"

echo "Date metadata saved to $date_file"
echo "Dates decimal metadata saved to $dates_decimal_file"
echo "Location metadata saved to $location_metadata"
echo "Compact location metadata saved to $location_file"