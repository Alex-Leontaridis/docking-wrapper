#!/usr/bin/env python3
"""
Test script to verify the new path manager and logging system.
"""

import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

def test_path_manager():
    """Test the path manager functionality."""
    print("Testing Path Manager...")
    
    try:
        from utils.path_manager import get_path_manager, get_path, get_absolute_path, ensure_dir
        
        # Get path manager
        path_manager = get_path_manager()
        print(f"✓ Path manager initialized")
        
        # Test basic path access
        project_root = get_path("project_root")
        print(f"✓ Project root: {project_root}")
        
        # Test directory creation
        test_dir = ensure_dir("temp", "test")
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
        return False

def test_logging():
    """Test the logging functionality."""
    print("\nTesting Logging System...")
    
    try:
        from utils.logging import setup_logging, log_startup, log_shutdown, log_error_with_context
        
        # Setup logger
        logger = setup_logging('TestLogger', level='INFO')
        print(f"✓ Logger initialized")
        
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
        return False

def test_integration():
    """Test integration between path manager and logging."""
    print("\nTesting Integration...")
    
    try:
        from utils.path_manager import get_path_manager
        from utils.logging import setup_logging
        
        # Setup with path manager integration
        path_manager = get_path_manager()
        log_dir = path_manager.ensure_dir("logs")
        
        if log_dir:
            log_file = log_dir / "test_integration.log"
            logger = setup_logging('IntegrationTest', str(log_file))
            
            logger.info("Integration test started")
            logger.info(f"Platform: {path_manager.get_platform_info()['system']}")
            logger.info(f"Available tools: {path_manager.list_available_tools()}")
            logger.info("Integration test completed")
            
            print(f"✓ Integration test passed - log file: {log_file}")
            return True
        else:
            print("✗ Could not create log directory")
            return False
            
    except Exception as e:
        print(f"✗ Integration test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("=" * 60)
    print("PATH MANAGER AND LOGGING SYSTEM TEST")
    print("=" * 60)
    
    # Run tests
    path_manager_ok = test_path_manager()
    logging_ok = test_logging()
    integration_ok = test_integration()
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    print(f"Path Manager: {'✓ PASSED' if path_manager_ok else '✗ FAILED'}")
    print(f"Logging System: {'✓ PASSED' if logging_ok else '✗ FAILED'}")
    print(f"Integration: {'✓ PASSED' if integration_ok else '✗ FAILED'}")
    
    if all([path_manager_ok, logging_ok, integration_ok]):
        print("\n🎉 ALL TESTS PASSED! The new system is working correctly.")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED. Please check the implementation.")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 