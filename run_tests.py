#!/usr/bin/env python3
"""
Test runner for the docking wrapper fixes.
"""

import unittest
import sys
import os
from pathlib import Path

def run_all_tests():
    """Run all tests and return results."""
    # Add tests directory to path
    tests_dir = Path(__file__).parent / 'tests'
    sys.path.insert(0, str(tests_dir))
    
    # Discover and run all tests
    loader = unittest.TestLoader()
    suite = loader.discover(str(tests_dir), pattern='test_*.py')
    
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result

def run_specific_test(test_name):
    """Run a specific test module."""
    tests_dir = Path(__file__).parent / 'tests'
    sys.path.insert(0, str(tests_dir))
    
    # Import and run specific test
    test_module = __import__(test_name)
    suite = unittest.TestLoader().loadTestsFromModule(test_module)
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result

def main():
    """Main test runner."""
    print("=" * 60)
    print("DOCKING WRAPPER FIXES - TEST SUITE")
    print("=" * 60)
    
    if len(sys.argv) > 1:
        # Run specific test
        test_name = sys.argv[1]
        print(f"Running specific test: {test_name}")
        result = run_specific_test(test_name)
    else:
        # Run all tests
        print("Running all tests...")
        result = run_all_tests()
    
    # Print summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped) if hasattr(result, 'skipped') else 0}")
    
    if result.failures:
        print("\nFAILURES:")
        for test, traceback in result.failures:
            print(f"  {test}: {traceback}")
    
    if result.errors:
        print("\nERRORS:")
        for test, traceback in result.errors:
            print(f"  {test}: {traceback}")
    
    # Return appropriate exit code
    if result.wasSuccessful():
        print("\n✅ All tests passed!")
        return 0
    else:
        print("\n❌ Some tests failed!")
        return 1

if __name__ == '__main__':
    sys.exit(main()) 