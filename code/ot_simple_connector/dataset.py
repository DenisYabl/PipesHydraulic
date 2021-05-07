# -*- coding: utf-8 -*-
import json
import logging


class Dataset:

    logger = logging.getLogger('ot_simple_connector')

    def __init__(self, session, schema_url, data_urls):
        self.session = session
        self.data_urls = data_urls
        self.schema_url = schema_url
        self._buffer = []
        self._point = 0
        self._schema = None

    @property
    def schema(self):
        if self._schema is None:
            self._schema = self.session.get(self.schema_url).text
        return self._schema

    def load(self):
        chunks = []
        for url in self.data_urls:
            chunk = self.session.get(url).text
            chunk = chunk.split('\n')
            chunk = [json.loads(chunk) for chunk in chunk if chunk]
            chunks += chunk
        return chunks

    def __iter__(self):
        return self

    def __next__(self):
        if self._buffer:
            return self._buffer.pop()
        elif self._point < len(self.data_urls):
            url = self.data_urls[self._point]
            self.logger.debug('Loading chunk %s.' % url)
            self._point += 1
            chunk = self.session.get(url)
            chunk = chunk.content
            chunk = chunk.decode('utf-8')
            self.logger.debug('Splitting chunk.')
            chunk = chunk.split('\n')
            self.logger.debug('Preparing chunk buffer.')
            chunk = [json.loads(chunk) for chunk in chunk if chunk]
            self._buffer = chunk
            self.logger.debug('Chunk buffer is ready.')
            return self.__next__()
        else:
            raise StopIteration

    def next(self):
        return self.__next__()
