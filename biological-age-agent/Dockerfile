# Dockerfile
FROM python:3.13-slim

WORKDIR /app

# 安装系统依赖（openpyxl 需要）
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# 安装 Python 依赖
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 复制代码（排除 .env、__pycache__ 等）
COPY . .

# 模型默认包含在镜像中
RUN test -f model/bio_age_prediction_model.pkl || echo "WARNING: model file not found"

# 暴露端口
EXPOSE 8000

# 启动
CMD ["python", "run.py"]
