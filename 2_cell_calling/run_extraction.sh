#!/bin/bash
# Helper script to show available CellBender conditions and run the extraction
# 
# Usage: 
#   ./run_extraction.sh                    # Interactive mode - shows menu
#   ./run_extraction.sh rmHs               # Run for rmHs condition
#   ./run_extraction.sh wHsPf              # Run for wHsPf condition
#   ./run_extraction.sh rmHs /path/to/out  # Specify output directory

BASE_PATH="/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/data/processed/Pf"
SCRIPT_DIR="/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/mali22_code_lnk/cell_calling"
DEFAULT_OUTPUT_DIR="${SCRIPT_DIR}"

# Parse output directory if provided
OUTPUT_DIR="${2:-$DEFAULT_OUTPUT_DIR}"

# Define conditions
declare -A CONDITIONS
CONDITIONS=(
    ["rmHs"]="${BASE_PATH}/MSC*/cbender_custom2_rmHs/cb_LR*/cb_report.html"
    ["wHsPf"]="${BASE_PATH}/MSC*/cbender_custom2_wHsPf/cb_LR*/cb_report.html"
)

# Function to count files for a pattern
count_files() {
    local pattern="$1"
    echo $(ls $pattern 2>/dev/null | wc -l)
}

# Function to show menu
show_menu() {
    echo "========================================"
    echo "CellBender Learning Curve Extraction"
    echo "========================================"
    echo ""
    echo "Available conditions:"
    echo ""
    
    local i=1
    for condition in "${!CONDITIONS[@]}"; do
        local pattern="${CONDITIONS[$condition]}"
        local count=$(count_files "$pattern")
        printf "%d) %-10s (%3d files)\n" $i "$condition" $count
        ((i++))
    done
    
    echo ""
    echo "0) Exit"
    echo ""
}

# Function to run extraction
run_extraction() {
    local condition="$1"
    local pattern="${CONDITIONS[$condition]}"
    
    if [ -z "$pattern" ]; then
        echo "Error: Unknown condition '$condition'"
        exit 1
    fi
    
    local count=$(count_files "$pattern")
    
    if [ "$count" -eq 0 ]; then
        echo "Warning: No files found for condition '$condition'"
        echo "Pattern: $pattern"
        exit 1
    fi
    
    echo "========================================"
    echo "Processing condition: $condition"
    echo "Files to process: $count"
    echo "Pattern: $pattern"
    echo "Output directory: $OUTPUT_DIR"
    echo "========================================"
    echo ""
    
    echo "Running R script..."
    Rscript "${SCRIPT_DIR}/extract_cellbender_learning_curve_assessment.R" "$pattern" "$OUTPUT_DIR"
}

# Main script
if [ $# -eq 0 ]; then
    # Interactive mode
    while true; do
        show_menu
        read -p "Enter choice: " choice
        
        case $choice in
            0)
                echo "Exiting..."
                exit 0
                ;;
            1)
                run_extraction "rmHs"
                break
                ;;
            2)
                run_extraction "wHsPf"
                break
                ;;
            *)
                echo "Invalid choice. Please try again."
                echo ""
                ;;
        esac
    done
else
    # Command line mode
    run_extraction "$1"
fi
