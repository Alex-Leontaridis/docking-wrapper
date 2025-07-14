import sentry_sdk

def init_sentry(dsn: str, environment: str = None):
    sentry_sdk.init(
        dsn=dsn,
        environment=environment,
        traces_sample_rate=1.0,  # Adjust as needed
        send_default_pii=True
    ) 