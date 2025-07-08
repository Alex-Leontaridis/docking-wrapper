#!/usr/bin/env python3
"""
Simple test script to verify the new path manager and logging system.
Tests modules separately to avoid circular import issues.
"""

import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

def test_path_manager_only():
    """Test only the path manager functionality."""
    print("Testing Path Manager (standalone)...")
    
    try:
        # Import path manager
        from utils.path_manager import PathManager
        
        # Create a new instance (not using global)
        path_manager = PathManager()
        print(f"✓ Path manager created")
        
        # Test basic path access
        project_root = path_manager.get_path("project_root")
        print(f"✓ Project root: {project_root}")
        
        # Test directory creation
        test_dir = path_manager.ensure_dir("temp", "test")
        if test_dir:
            print(f"✓ Test directory created: {test_dir}")
        
        # Test platform info
        platform_info = path_manager.get_platform_info()
        print(f"✓ Platform: {platform_info['system']} {platform_info['release']}")
        
        # Test tool availability
        tools = path_manager.list_available_tools()
        print(f"✓ Available tools: {tools}")
        
        print("✓ Path manager tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Path manager test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_logging_only():
    """Test only the logging functionality."""
    print("\nTesting Logging System (standalone)...")
    
    try:
        # Import logging
        from utils.logging import DockingLogger
        
        # Create logger without path manager integration
        logger = DockingLogger('TestLogger', level='INFO')
        print(f"✓ Logger created")
        
        # Test basic logging
        logger.info("Test info message")
        logger.warning("Test warning message")
        logger.error("Test error message", error_code='FILE_NOT_FOUND')
        print(f"✓ Basic logging works")
        
        # Test performance tracking
        with logger.timer("test_operation"):
            import time
            time.sleep(0.1)  # Simulate work
        print(f"✓ Performance tracking works")
        
        # Test convenience functions
        from utils.logging import log_startup, log_shutdown, log_error_with_context
        
        log_startup('TestModule', '1.0.0')
        log_shutdown('TestModule')
        print(f"✓ Convenience functions work")
        
        # Test error context
        try:
            raise ValueError("Test error")
        except Exception as e:
            log_error_with_context(e, "During testing", error_code='SYSTEM_ERROR')
        print(f"✓ Error context logging works")
        
        print("✓ Logging tests passed")
        return True
        
    except Exception as e:
        print(f"✗ Logging test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_basic_integration():
    """Test basic integration without circular imports."""
    print("\nTesting Basic Integration...")
    
    try:
        # Test that both modules can be imported
        from utils import path_manager
        from utils import logging
        
        print(f"✓ Both modules imported successfully")
        
        # Test that path manager can be used
        pm = path_manager.PathManager()
        print(f"✓ Path manager instance created")
        
        # Test that logger can be used
        logger = logging.DockingLogger('IntegrationTest')
        print(f"✓ Logger instance created")
        
        # Test basic functionality
        project_root = pm.get_path("project_root")
        logger.info(f"Project root: {project_root}")
        
        print(f"✓ Basic integration test passed")
        return True
        
    except Exception as e:
        print(f"✗ Basic integration test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("SIMPLE PATH MANAGER AND LOGGING SYSTEM TEST")
    print("=" * 60)
    
    # Run tests
    path_manager_ok = test_path_manager_only()
    logging_ok = test_logging_only()
    integration_ok = test_basic_integration()
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    print(f"Path Manager: {'✓ PASSED' if path_manager_ok else '✗ FAILED'}")
    print(f"Logging System: {'✓ PASSED' if logging_ok else '✗ FAILED'}")
    print(f"Basic Integration: {'✓ PASSED' if integration_ok else '✗ FAILED'}")
    
    if all([path_manager_ok, logging_ok, integration_ok]):
        print("\n🎉 ALL TESTS PASSED! The new system is working correctly.")
        print("\nNext steps:")
        print("1. The path manager and logging systems are working independently")
        print("2. You can now use them in your scripts")
        print("3. Run the batch pipeline to test the full integration")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED. Please check the implementation.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 