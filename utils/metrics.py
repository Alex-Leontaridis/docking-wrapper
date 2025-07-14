from prometheus_client import Histogram, Counter, start_http_server

# Histograms for durations
DOCKING_DURATION = Histogram('docking_duration_seconds', 'Time spent in docking', ['engine'])
INFERENCE_DURATION = Histogram('inference_duration_seconds', 'Time spent in inference', ['engine'])

# Counter for errors
PIPELINE_ERRORS = Counter('pipeline_errors_total', 'Total pipeline errors', ['stage'])

# Start the Prometheus metrics server (default port 8000)
def start_metrics_server(port=8000):
    start_http_server(port) 