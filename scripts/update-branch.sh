#!/bin/bash

# Update Branch Script
# Usage: ./scripts/update-branch.sh [branch-name]
# If no branch name provided, updates current branch

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

# Check if we're in a git repository
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    print_error "Not in a git repository. Please run this script from the repository root."
    exit 1
fi

# Get current branch
CURRENT_BRANCH=$(git branch --show-current)

# Determine target branch
if [ $# -eq 1 ]; then
    TARGET_BRANCH=$1
else
    TARGET_BRANCH=$CURRENT_BRANCH
fi

print_status "Updating branch: $TARGET_BRANCH"

# Check if target branch exists
if ! git show-ref --verify --quiet refs/heads/$TARGET_BRANCH; then
    print_error "Branch $TARGET_BRANCH doesn't exist locally"
    exit 1
fi

# Check if we have uncommitted changes
if ! git diff-index --quiet HEAD --; then
    print_warning "You have uncommitted changes. Stashing them..."
    git stash push -m "Auto-stash before updating branch"
    STASHED=true
fi

# Fetch latest changes
print_status "Fetching latest changes from remote..."
git fetch origin

# Check if develop branch exists
if ! git show-ref --verify --quiet refs/heads/develop; then
    print_error "Develop branch doesn't exist. Please create it first."
    exit 1
fi

# Update develop branch
print_status "Updating develop branch..."
git checkout develop
git pull origin develop

# Switch to target branch
print_status "Switching to target branch..."
git checkout $TARGET_BRANCH

# Rebase on develop
print_status "Rebasing on develop..."
if git rebase develop; then
    print_success "Successfully updated $TARGET_BRANCH with latest changes from develop"
else
    print_error "Rebase failed. You may need to resolve conflicts manually."
    print_status "After resolving conflicts, run: git rebase --continue"
    exit 1
fi

# Push changes if branch exists on remote
if git ls-remote --heads origin $TARGET_BRANCH | grep -q $TARGET_BRANCH; then
    print_status "Pushing updated branch to remote..."
    if git push origin $TARGET_BRANCH --force-with-lease; then
        print_success "Successfully pushed updated branch to remote"
    else
        print_warning "Failed to push to remote. You may need to force push manually."
        print_status "Run: git push origin $TARGET_BRANCH --force-with-lease"
    fi
fi

# Restore stashed changes if any
if [ "$STASHED" = true ]; then
    print_status "Restoring stashed changes..."
    git stash pop
fi

print_success "Branch update completed successfully!" 