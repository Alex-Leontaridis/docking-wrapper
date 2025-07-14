#!/bin/bash

# Create Feature Branch Script
# Usage: ./scripts/create-feature.sh <component> <description>
# Example: ./scripts/create-feature.sh equibind integration

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

# Check arguments
if [ $# -lt 2 ]; then
    print_error "Usage: $0 <component> <description>"
    print_error "Example: $0 equibind integration"
    print_error "Example: $0 gnina enhancement"
    exit 1
fi

COMPONENT=$1
shift
DESCRIPTION="$*"

# Convert to lowercase and replace spaces with hyphens
COMPONENT=$(echo "$COMPONENT" | tr '[:upper:]' '[:lower:]' | tr ' ' '-')
DESCRIPTION=$(echo "$DESCRIPTION" | tr '[:upper:]' '[:lower:]' | tr ' ' '-')

BRANCH_NAME="feature/${COMPONENT}-${DESCRIPTION}"

print_status "Creating feature branch: $BRANCH_NAME"

# Check if branch already exists
if git show-ref --verify --quiet refs/heads/$BRANCH_NAME; then
    print_error "Branch $BRANCH_NAME already exists locally"
    exit 1
fi

if git ls-remote --heads origin $BRANCH_NAME | grep -q $BRANCH_NAME; then
    print_error "Branch $BRANCH_NAME already exists on remote"
    exit 1
fi

# Fetch latest changes
print_status "Fetching latest changes from remote..."
git fetch origin

# Check if develop branch exists
if ! git show-ref --verify --quiet refs/heads/develop; then
    print_warning "Develop branch doesn't exist. Creating from main..."
    git checkout main
    git pull origin main
    git checkout -b develop
    git push -u origin develop
fi

# Switch to develop and update
print_status "Updating develop branch..."
git checkout develop
git pull origin develop

# Create feature branch
print_status "Creating feature branch from develop..."
git checkout -b $BRANCH_NAME

# Push to remote
print_status "Pushing feature branch to remote..."
git push -u origin $BRANCH_NAME

print_success "Feature branch '$BRANCH_NAME' created successfully!"
print_status "You can now start developing your feature."
print_status "When ready, create a Pull Request from $BRANCH_NAME to develop"

# Show next steps
echo ""
print_status "Next steps:"
echo "  1. Start developing your feature"
echo "  2. Make regular commits with descriptive messages"
echo "  3. Push changes: git push origin $BRANCH_NAME"
echo "  4. Create Pull Request to develop when ready"
echo "  5. Clean up branch after merge: git branch -d $BRANCH_NAME" 