FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

# 容器内服务监听端口（与 --port 保持一致）
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD python -c "import urllib.request,sys; sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:8000/').status == 200 else 1)"

# 端口可通过 uvicorn 参数覆盖，例如：
#   docker run -p 9000:9000 <image> uvicorn app.main:app --host 0.0.0.0 --port 9000
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
