# -*- coding: utf-8 -*-
import logging
import time
import uuid
from ot_simple_connector.job import Job


class Jobs:

    logger = logging.getLogger('ot_simple_connector')

    def __init__(self, session):
        self.session = session

    def create(self, query_text, cache_ttl=0, tws=0, twf=0, sid=uuid.uuid1(), blocking=True, timeout=0):
        job = Job(self.session, query_text, cache_ttl, tws, twf, sid)
        job.create()
        if blocking:
            job_time_start = time.time()
            job_status = job.status
            self.logger.debug('Job[%s] is %s.' % (job.sid, job_status))
            while job_status not in ['failed', 'success']:
                if timeout and time.time() - job_time_start > timeout:
                    raise Exception('Job failed because of timeout.')
                time.sleep(1)
                job_status = job.status
                self.logger.debug('Job[%s] is %s.' % (job.sid, job_status))
                if job_status == 'failed':
                    raise Exception('Job[%s] failed because of %s.' % (job.sid, job.msg))
        return job
