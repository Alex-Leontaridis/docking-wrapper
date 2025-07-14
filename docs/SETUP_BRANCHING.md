# Setting Up the Branching Strategy

This guide provides step-by-step instructions for implementing the structured branching model in your repository.

## Prerequisites

- Git installed and configured
- Access to the repository with admin privileges
- Basic understanding of Git workflows

## Step 1: Initial Repository Setup

### 1.1 Configure Git Hooks

Make the commit message hook executable:
```bash
chmod +x .git/hooks/commit-msg
```

### 1.2 Set Commit Message Template

Configure Git to use the commit message template:
```bash
git config commit.template .gitmessage
```

### 1.3 Configure Git Aliases (Optional)

Add useful Git aliases to your `.gitconfig`:
```bash
git config --global alias.st status
git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.ci commit
git config --global alias.lg "log --oneline --graph --all"
git config --global alias.unstage "reset HEAD --"
git config --global alias.last "log -1 HEAD"
```

## Step 2: Create Develop Branch

### 2.1 Create and Push Develop Branch

```bash
# Ensure you're on main and it's up to date
git checkout main
git pull origin main

# Create develop branch
git checkout -b develop

# Push develop branch to remote
git push -u origin develop
```

### 2.2 Set Up Branch Protection (GitHub/GitLab)

#### GitHub Settings:
1. Go to repository Settings → Branches
2. Add rule for `main` branch:
   - Require pull request reviews before merging
   - Require status checks to pass before merging
   - Require branches to be up to date before merging
   - Restrict pushes that create files larger than 100 MB
   - Require linear history
3. Add rule for `develop` branch:
   - Require pull request reviews before merging
   - Require status checks to pass before merging
   - Allow force pushes for administrators only

#### GitLab Settings:
1. Go to repository Settings → Repository → Protected Branches
2. Protect `main` branch:
   - Allowed to merge: Maintainers
   - Allowed to push: No one
3. Protect `develop` branch:
   - Allowed to merge: Developers
   - Allowed to push: Maintainers

## Step 3: Set Up CI/CD Pipeline

### 3.1 Create GitHub Actions Workflow

Create `.github/workflows/ci.yml`:

```yaml
name: CI/CD Pipeline

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]

jobs:
  test:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.9'
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
    
    - name: Run tests
      run: |
        python -m pytest tests/
    
    - name: Run linting
      run: |
        pip install flake8 black
        flake8 scripts/ --max-line-length=88
        black --check scripts/

  security:
    runs-on: ubuntu-latest
    
    steps:
    - uses: actions/checkout@v3
    
    - name: Run security scan
      uses: github/codeql-action/init@v2
      with:
        languages: python
    
    - name: Perform CodeQL Analysis
      uses: github/codeql-action/analyze@v2
```

### 3.2 Create GitLab CI Pipeline

Create `.gitlab-ci.yml`:

```yaml
stages:
  - test
  - security
  - deploy

test:
  stage: test
  image: python:3.9
  script:
    - pip install -r requirements.txt
    - python -m pytest tests/
    - pip install flake8 black
    - flake8 scripts/ --max-line-length=88
    - black --check scripts/
  only:
    - main
    - develop
    - merge_requests

security:
  stage: security
  image: python:3.9
  script:
    - pip install bandit
    - bandit -r scripts/ -f json -o security-report.json
  artifacts:
    reports:
      security: security-report.json
  only:
    - main
    - develop
    - merge_requests

deploy_staging:
  stage: deploy
  script:
    - echo "Deploying to staging"
  environment:
    name: staging
  only:
    - develop

deploy_production:
  stage: deploy
  script:
    - echo "Deploying to production"
  environment:
    name: production
  only:
    - main
```

## Step 4: Create Initial Feature Branches

### 4.1 Create Example Feature Branches

```bash
# EquiBind integration
./scripts/create-feature.sh equibind integration

# GNINA enhancement
./scripts/create-feature.sh gnina enhancement

# Docker containerization
./scripts/create-feature.sh docker containerization

# Testing framework
./scripts/create-feature.sh testing framework
```

### 4.2 Set Up Development Environment

Create a development environment setup script:

```bash
# scripts/setup-dev.sh
#!/bin/bash

echo "Setting up development environment..."

# Install dependencies
pip install -r requirements.txt

# Install development dependencies
pip install pytest pytest-cov flake8 black mypy

# Set up pre-commit hooks
pip install pre-commit
pre-commit install

# Create necessary directories
mkdir -p logs outputs test_outputs

echo "Development environment setup complete!"
```

## Step 5: Documentation Setup

### 5.1 Update README

Add a section to the main README about the branching strategy:

```markdown
## Development Workflow

This project follows a structured branching model:

- `main`: Production-ready code
- `develop`: Integration branch for feature development
- `feature/*`: Short-lived branches for new features
- `release/*`: Release preparation branches
- `hotfix/*`: Critical bug fixes

For detailed information, see:
- [Branching Strategy](docs/BRANCHING_STRATEGY.md)
- [Developer Workflow](docs/DEVELOPER_WORKFLOW.md)
- [Setup Guide](docs/SETUP_BRANCHING.md)

### Quick Start for Developers

1. Create a feature branch:
   ```bash
   ./scripts/create-feature.sh <component> <description>
   ```

2. Develop your feature and commit changes:
   ```bash
   git add .
   git commit -m "feat(component): add new functionality"
   ```

3. Push and create a Pull Request:
   ```bash
   git push origin feature/your-branch
   ```

4. After review and approval, merge to `develop`
```

### 5.2 Create Team Guidelines

Create `docs/TEAM_GUIDELINES.md`:

```markdown
# Team Development Guidelines

## Code Review Process

1. All changes must go through Pull Request review
2. At least one approval required for feature branches
3. At least two approvals required for main/develop branches
4. All CI checks must pass before merging

## Commit Message Standards

- Use conventional commit format
- Keep messages under 72 characters
- Reference issues when applicable
- Use imperative mood

## Branch Naming

- Feature branches: `feature/component-description`
- Release branches: `release/v1.0.0`
- Hotfix branches: `hotfix/issue-description`

## Testing Requirements

- Unit tests for new functionality
- Integration tests for API changes
- Manual testing for UI changes
- Performance testing for optimization changes
```

## Step 6: Automation Setup

### 6.1 Create Automated Scripts

#### Branch Cleanup Automation

Create a GitHub Action for automatic branch cleanup:

```yaml
# .github/workflows/cleanup-branches.yml
name: Cleanup Merged Branches

on:
  schedule:
    - cron: '0 2 * * 0'  # Weekly on Sunday at 2 AM

jobs:
  cleanup:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
      with:
        fetch-depth: 0
    
    - name: Cleanup merged branches
      run: |
        git config --global user.email "action@github.com"
        git config --global user.name "GitHub Action"
        
        # Delete merged feature branches
        git branch --merged origin/develop | grep 'feature/' | xargs -n 1 git push origin --delete || true
```

#### Release Automation

Create a release automation script:

```bash
# scripts/create-release.sh
#!/bin/bash

VERSION=$1

if [ -z "$VERSION" ]; then
    echo "Usage: $0 <version>"
    echo "Example: $0 1.0.0"
    exit 1
fi

# Create release branch
git checkout develop
git pull origin develop
git checkout -b release/v$VERSION

# Update version in files
echo $VERSION > VERSION
sed -i "s/version = .*/version = \"$VERSION\"/" setup.py

# Commit version bump
git add VERSION setup.py
git commit -m "chore(release): bump version to $VERSION"

# Push release branch
git push origin release/v$VERSION

echo "Release branch created: release/v$VERSION"
echo "Create Pull Request from release/v$VERSION to main"
```

## Step 7: Training and Adoption

### 7.1 Team Training

1. **Workshop Session**: Conduct a 2-hour workshop covering:
   - Branching strategy overview
   - Hands-on practice with scripts
   - Common scenarios and troubleshooting

2. **Documentation Review**: Ensure all team members have read:
   - Branching Strategy document
   - Developer Workflow guide
   - Team Guidelines

3. **Practice Exercises**: Create practice scenarios:
   - Creating feature branches
   - Resolving merge conflicts
   - Creating releases
   - Handling hotfixes

### 7.2 Migration Checklist

- [ ] Develop branch created and protected
- [ ] CI/CD pipeline configured
- [ ] Branch protection rules set
- [ ] Documentation updated
- [ ] Team trained on new workflow
- [ ] Scripts tested and working
- [ ] Example feature branches created
- [ ] Release process tested

### 7.3 Monitoring and Maintenance

#### Regular Tasks:
- Weekly branch cleanup
- Monthly workflow review
- Quarterly strategy assessment

#### Metrics to Track:
- Time from feature branch to merge
- Number of merge conflicts
- Release frequency
- Code review turnaround time

## Troubleshooting

### Common Issues

#### Branch Protection Conflicts
```bash
# If you need to bypass protection temporarily
git push origin branch-name --force-with-lease
```

#### Merge Conflicts
```bash
# During rebase
git add <resolved-files>
git rebase --continue

# During merge
git add <resolved-files>
git commit
```

#### Hook Issues
```bash
# Skip hooks if needed (not recommended)
git commit --no-verify -m "message"

# Debug hook issues
chmod +x .git/hooks/commit-msg
```

### Getting Help

1. Check the documentation first
2. Review the troubleshooting section
3. Ask in team chat
4. Create an issue for persistent problems
5. Contact repository maintainers

## Next Steps

After implementing the branching strategy:

1. **Monitor Adoption**: Track how well the team is following the new workflow
2. **Gather Feedback**: Collect feedback from team members
3. **Iterate**: Make improvements based on feedback
4. **Automate**: Add more automation as needed
5. **Scale**: Adapt the strategy as the team grows 