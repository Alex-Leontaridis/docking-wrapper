#!/bin/bash

# Cleanup Branches Script
# Usage: ./scripts/cleanup-branches.sh [--dry-run] [--remote]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse arguments
DRY_RUN=false
CLEANUP_REMOTE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --remote)
            CLEANUP_REMOTE=true
            shift
            ;;
        *)
            print_error "Unknown option: $1"
            print_error "Usage: $0 [--dry-run] [--remote]"
            exit 1
            ;;
    esac
done

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    print_error "Not in a git repository. Please run this script from the repository root."
    exit 1
fi

print_status "Starting branch cleanup..."

# Fetch latest changes
print_status "Fetching latest changes from remote..."
git fetch origin --prune

# Get current branch
CURRENT_BRANCH=$(git branch --show-current)

# Function to cleanup local branches
cleanup_local_branches() {
    print_status "Cleaning up local branches..."
    
    # Get list of merged branches
    MERGED_BRANCHES=$(git branch --merged develop | grep -v "develop" | grep -v "main" | grep -v "master" | sed 's/^[[:space:]]*//')
    
    if [ -z "$MERGED_BRANCHES" ]; then
        print_status "No merged branches to clean up locally."
        return
    fi
    
    echo "Merged branches to be deleted:"
    echo "$MERGED_BRANCHES"
    echo ""
    
    if [ "$DRY_RUN" = true ]; then
        print_warning "DRY RUN: Would delete the following local branches:"
        echo "$MERGED_BRANCHES"
        return
    fi
    
    # Confirm deletion
    read -p "Do you want to delete these merged branches? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_status "Branch cleanup cancelled."
        return
    fi
    
    # Delete merged branches
    for branch in $MERGED_BRANCHES; do
        if [ "$branch" != "$CURRENT_BRANCH" ]; then
            if git branch -d "$branch" 2>/dev/null; then
                print_success "Deleted local branch: $branch"
            else
                print_warning "Could not delete local branch: $branch (may have unmerged changes)"
            fi
        else
            print_warning "Skipping current branch: $branch"
        fi
    done
}

# Function to cleanup remote branches
cleanup_remote_branches() {
    if [ "$CLEANUP_REMOTE" != true ]; then
        return
    fi
    
    print_status "Cleaning up remote branches..."
    
    # Get list of remote branches that have been merged
    REMOTE_MERGED_BRANCHES=$(git branch -r --merged origin/develop | grep -v "origin/develop" | grep -v "origin/main" | grep -v "origin/master" | sed 's/origin\///')
    
    if [ -z "$REMOTE_MERGED_BRANCHES" ]; then
        print_status "No merged remote branches to clean up."
        return
    fi
    
    echo "Merged remote branches to be deleted:"
    echo "$REMOTE_MERGED_BRANCHES"
    echo ""
    
    if [ "$DRY_RUN" = true ]; then
        print_warning "DRY RUN: Would delete the following remote branches:"
        echo "$REMOTE_MERGED_BRANCHES"
        return
    fi
    
    # Confirm deletion
    read -p "Do you want to delete these merged remote branches? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_status "Remote branch cleanup cancelled."
        return
    fi
    
    # Delete remote branches
    for branch in $REMOTE_MERGED_BRANCHES; do
        if git push origin --delete "$branch" 2>/dev/null; then
            print_success "Deleted remote branch: $branch"
        else
            print_warning "Could not delete remote branch: $branch"
        fi
    done
}

# Function to show stale branches
show_stale_branches() {
    print_status "Checking for stale branches..."
    
    # Show branches that haven't been updated in 30 days
    STALE_BRANCHES=$(git for-each-ref --format='%(refname:short) %(committerdate:relative)' refs/heads | grep -v "develop" | grep -v "main" | grep -v "master" | grep -E "(days|weeks|months) ago" | head -10)
    
    if [ -n "$STALE_BRANCHES" ]; then
        print_warning "Potentially stale branches (not updated recently):"
        echo "$STALE_BRANCHES"
        echo ""
    fi
}

# Main cleanup process
cleanup_local_branches
cleanup_remote_branches
show_stale_branches

print_success "Branch cleanup completed!"

if [ "$DRY_RUN" = true ]; then
    print_warning "This was a dry run. No branches were actually deleted."
fi 