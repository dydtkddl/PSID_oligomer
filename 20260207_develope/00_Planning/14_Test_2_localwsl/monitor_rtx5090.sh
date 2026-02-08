#!/bin/bash
# monitor_rtx5090.sh - RTX 5090 전용 모니터링

trap "tput cnorm; clear; exit" SIGINT SIGTERM
tput civis

while true; do
    clear
    echo "=========================================="
    echo "🚀 RTX 5090 실시간 모니터링"
    echo "시간: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    
    # GPU 상세 정보
    nvidia-smi --query-gpu=name,temperature.gpu,power.draw,power.limit,utilization.gpu,memory.used,memory.total \
        --format=csv,noheader | \
        awk -F', ' '{
            printf "GPU 모델: %s\n", $1
            printf "온도    : %s°C\n", $2
            printf "전력    : %s / %s W\n", $3, $4
            printf "사용률  : %s%%\n", $5
            printf "메모리  : %s / %s\n", $6, $7
            
            # 메모리 사용률 바 표시
            split($6, used, " "); split($7, total, " ");
            if (total[1] > 0) {
                pct = (used[1] / total[1]) * 100;
                printf "메모리바: [";
                for(i=0;i<50;i++) {
                    if(i < pct*50/100) printf "█"; else printf "░";
                }
                printf "] %.1f%%\n", pct;
            }
        }'
    
    echo "=========================================="
    echo "💡 팁: RTX 5090은 32GB VRAM으로 초대규모 시뮬레이션 가능"
    echo "Ctrl+C로 종료"
    
    sleep 2
done
