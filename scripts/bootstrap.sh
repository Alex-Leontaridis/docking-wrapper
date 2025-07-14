#!/bin/bash

# Set environment variables for observability
export SENTRY_DSN=${SENTRY_DSN:-""}
export METRICS_PORT=${METRICS_PORT:-8000}

# Start Prometheus metrics server in the background
python -c "from utils.metrics import start_metrics_server; start_metrics_server(port=int(\"${METRICS_PORT}\"))" & 