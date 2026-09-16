FROM python:3.11-slim
ENV MPLBACKEND=Agg PYTHONHASHSEED=17 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1
WORKDIR /app
COPY pyproject.toml README.md requirements-tested.txt ./
COPY scrna_workflow ./scrna_workflow
RUN pip install --no-cache-dir -r requirements-tested.txt && pip install --no-cache-dir --no-deps .
ENTRYPOINT ["python", "-m", "scrna_workflow"]
