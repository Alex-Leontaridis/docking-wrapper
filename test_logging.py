from utils.logging import setup_logging, ERROR_CODES
logger = setup_logging('TestLogger', 'test.log', 'DEBUG')
logger.info('Test info message')
logger.error('Test error message', error_code='FILE_NOT_FOUND')
logger.warning('Test warning message')
print('Available error codes:')
for code, number in ERROR_CODES.items():
    print(f'  {code}: {number}') 